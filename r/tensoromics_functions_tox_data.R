dyn.load("build/libtensor-omics.so")
source("r/tensoromics_functions.R")
source("r/error_handling.R")

debug <- FALSE

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

#' Read expression values from multiple tabular (csv/tsv) files into a numeric matrix
#' @param file_list Character vector of filenames
#' @param gene_ids Character vector of gene IDs to match
#' @param n_header_rows Number of header rows to skip in each file
#' @param gene_col Column index (1-based) of gene IDs in each file
#' @param value_cols Vector of column indices (1-based) of expression values in each file
#' @param start_row Row index (1-based) to start reading data in each file
#' @param delimiter Delimiter used in the files (default: tab)
#' @param n_samples Total number of samples (rows) to read across all files
#' @return A list with elements:
#'   - expression_vectors: A numeric matrix of dimensions (n_samples x n_genes)
#'   - ierr: Integer error code (0 if successful)
read_expression_vectors_tsv <- function(file_list, gene_ids, 
                                    n_header_rows, gene_col, value_cols, delimiter = "\t", n_samples) {
    nfiles <- length(file_list)
    ngenes <- length(gene_ids)
    gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator

    # Convert file names to raw bytes
    file_raw_list <- lapply(file_list, function(f) {
      raw_bytes <- charToRaw(enc2native(f))
      # Pad with zeros to fixed length of 256
      if (length(raw_bytes) < 256) {
        raw_bytes <- c(raw_bytes, raw(256 - length(raw_bytes)))
      }
      raw_bytes[1:256]
    })
    file_raw <- do.call(cbind, file_raw_list)

    # Convert gene IDs to raw matrix
    gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
    
    # Convert delimiter to raw
    delimiter_raw <- charToRaw(delimiter)
    
    # Prepare output
    expression_vectors <- double(n_samples * ngenes)
    
    out <- .Fortran("read_expression_vectors_tsv_R",
        file_list_raw = file_raw,              # Pass raw bytes directly
        file_list_len = as.integer(256),       # Fixed length of 256
        n_files = as.integer(nfiles),
        gene_ids_raw = gene_raw,               # Pass raw bytes directly
        gene_ids_len = as.integer(gene_len),
        n_genes = as.integer(ngenes),
        expression_vectors = as.double(expression_vectors),
        n_samples = as.integer(n_samples),
        n_header_rows = as.integer(n_header_rows),
        gene_col = as.integer(gene_col),
        value_cols = as.integer(value_cols),
        n_value_cols = as.integer(length(value_cols)),
        ierr = 0L,
        delimiter_raw = delimiter_raw,         # Pass raw bytes directly
        dlen = as.integer(length(delimiter_raw))
    )

    check_err_code(out$ierr)
    
    list(
        expression_vectors = matrix(out$expression_vectors, nrow = n_samples, ncol = ngenes),
        ierr = out$ierr
    )
}

#' Read gene IDs from a single file
#' @param filename Character string of the filename
#' @param ngenes Number of gene IDs to read
#' @param gene_len Maximum length of each gene ID
#' @param n_header_rows Number of header rows to skip
#' @param gene_col Column index (1-based) of gene IDs in the file
#' @return A list with elements:
#'   - gene_ids: Character vector of gene IDs read from the file
#'   - ierr: Integer error code (0 if successful)
read_gene_ids_from_tsv_file <- function(filename, ngenes, gene_len, n_header_rows = 1, gene_col = 1) {
  # Convert filename to raw bytes
  filename_raw <- charToRaw(filename)
  
  gene_raw <- matrix(raw(1), nrow = gene_len + 1, ncol = ngenes)  # +1 for null terminator
  ierr <- integer(1)
  
  out <- .Fortran("read_gene_ids_from_tsv_file_R",
    filename_raw = filename_raw,              # Pass raw bytes directly
    fn_len = as.integer(length(filename_raw)),
    gene_ids_raw = gene_raw,                  # Pass raw bytes directly
    gene_ids_len = as.integer(gene_len + 1),  # +1 for null terminator
    n_genes = as.integer(ngenes),
    n_header_rows = as.integer(n_header_rows),
    gene_col = as.integer(gene_col),
    ierr = 0
  )
  check_err_code(out$ierr)
  
  list(
    gene_ids = raw_matrix_to_strings(out$gene_ids_raw, gene_len + 1),
    ierr = out$ierr
  )
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
  ngenes <- length(gene_ids)
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  
  # Convert filename to raw bytes
  filename_raw <- charToRaw(filename)
  
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  family_raw <- matrix(raw(1), nrow = family_len + 1, ncol = n_families)  # +1 for null terminator
  gene_to_fam <- integer(ngenes)
  ierr <- integer(1)
  
  out <- .Fortran("read_orthofinder_file_R",
    filename_raw = filename_raw,              # Pass raw bytes directly
    fn_len = as.integer(length(filename_raw)),
    gene_ids_raw = gene_raw,                  # Pass raw bytes directly
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    family_ids_raw = family_raw,              # Pass raw bytes directly
    family_ids_len = as.integer(family_len + 1),  # +1 for null terminator
    n_families = as.integer(n_families),
    gene_to_fam = as.integer(gene_to_fam),
    ierr = 0
  )
  check_err_code(out$ierr)
  list(
    family_ids = raw_matrix_to_strings(matrix(out$family_ids_raw, nrow = family_len + 1, ncol = n_families)),
    gene_to_fam = out$gene_to_fam,
    ierr = out$ierr
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
  # Erzeuge logische Maske: TRUE = behalten, FALSE = verwerfen
  mask <- gene_to_fam != 0L

  if(length(gene_ids) != ncol(expression_vectors) || length(gene_ids) != length(gene_to_fam)) {
    stop("Dimension mismatch: gene_ids, expression_vectors, and gene_to_fam must have compatible lengths")
  }
  
  # Filtere alle Arrays mit der gleichen Maske
  gene_ids <- gene_ids[mask]
  gene_to_fam <- gene_to_fam[mask]
  expression_vectors <- expression_vectors[, mask, drop = FALSE]
  
  n_genes_kept <- sum(mask)
  
  list(
    gene_ids = gene_ids,
    expression_vectors = expression_vectors,
    gene_to_fam = gene_to_fam,
    n_genes_kept = n_genes_kept
  )
}


# R wrappers for Fortran validation routines
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
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  fam_len  <- max(nchar(gene_family_ids)) + 1  # +1 for null terminator
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  fam_raw  <- strings_to_raw_matrix(gene_family_ids, fam_len)
  ierr <- integer(1)

  out <- .Fortran("validate_data_structure_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  gene_raw,                    # Pass raw bytes directly
                  as.integer(gene_len),
                  fam_raw,                     # Pass raw bytes directly
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Validate gene to family mapping
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param n_genes Number of genes
#' @param n_families Number of gene families
validate_gene_to_family_mapping <- function(gene_to_fam, n_genes, n_families) {
  ierr <- integer(1)
  out <- .Fortran("validate_gene_to_family_mapping_R",
                  as.integer(gene_to_fam),
                  as.integer(n_genes),
                  as.integer(n_families),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Validate expression data
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param check_non_negative Logical flag to check for non-negative values
validate_expression_data <- function(expression_vectors, n_genes, n_samples, check_non_negative = TRUE) {
  ierr <- integer(1)
  out <- .Fortran("validate_expression_data_R",
                  as.double(expression_vectors),
                  as.integer(n_genes),
                  as.integer(n_samples),
                  as.logical(check_non_negative),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Validate family centroids
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param n_families Number of gene families
#' @param n_samples Number of samples
validate_family_centroids <- function(family_centroids, n_families, n_samples) {
  ierr <- integer(1)
  out <- .Fortran("validate_family_centroids_R",
                  as.double(family_centroids),
                  as.integer(n_families),
                  as.integer(n_samples),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Validate shift vectors
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param n_samples Number of samples
validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples) {
  ierr <- integer(1)
  
  # Get dimensions
  n_genes <- ncol(expression_vectors)
  n_families <- ncol(family_centroids)
  
  out <- .Fortran("validate_shift_vectors_R",
                  as.double(shift_vectors),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.integer(gene_to_fam),
                  as.integer(n_genes),
                  as.integer(n_samples),
                  as.integer(n_families),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Validate uniqueness of string array
#' @param string_arr Character vector of gene IDs
#' @param n_strings Number of genes
#' Note: Uses hashset internally which may increase memory usage temporarily for large datasets
validate_string_array_uniqueness <- function(string_arr, n_strings) {
  str_len <- max(nchar(string_arr)) + 1  # +1 for null terminator
  str_raw <- strings_to_raw_matrix(string_arr, str_len)
  ierr <- integer(1)
  out <- .Fortran("validate_string_array_uniqueness_R",
                  str_raw,                    # Pass raw bytes directly
                  as.integer(str_len),
                  as.integer(n_strings),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
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
  # Determine max lengths for string encoding (+1 for null terminator)
  gene_len <- max(nchar(gene_ids)) + 1
  fam_len  <- max(nchar(gene_family_ids)) + 1

  # Convert to raw matrices
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  fam_raw  <- strings_to_raw_matrix(gene_family_ids, fam_len)

  ierr <- integer(1)

  out <- .Fortran("validate_all_data_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  gene_raw,                    # Pass raw bytes directly
                  as.integer(gene_len),
                  fam_raw,                     # Pass raw bytes directly
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)

  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

#' Low-level function to create zip archive from keys and filenames.
#' Directly calls the Fortran function.
#'
#' @param zip_filename Name of the zip file to create
#' @param keys Vector of keys for the manifest
#' @param filenames Vector of filenames to include in archive
create_zip_archive <- function(zip_filename, keys, filenames) {
  if (length(keys) != length(filenames)) {
    stop("Keys and filenames must have the same length")
  }
  
  # Convert to raw arrays for .Fortran call
  zip_raw <- charToRaw(zip_filename)
  zip_len <- length(zip_raw)
  
  # Calculate maximum lengths including null terminators
  max_key_len <- max(nchar(keys, type = "bytes")) + 1  # +1 for null terminator
  max_filename_len <- max(nchar(filenames, type = "bytes")) + 1
  count <- length(keys)
  
  # Create matrices initialized with null bytes
  keys_array <- matrix(as.raw(0), nrow = max_key_len, ncol = count)
  filenames_array <- matrix(as.raw(0), nrow = max_filename_len, ncol = count)
  
  # Fill arrays with proper null termination
  for (i in 1:count) {
    # Convert to raw and ensure null termination
    key_raw <- charToRaw(substr(keys[i], 1, max_key_len))
    filename_raw <- charToRaw(substr(filenames[i], 1, max_filename_len))
    
    # Check bounds and copy
    if (length(key_raw) > max_key_len) {
      stop("Key '", keys[i], "' is too long. Maximum length is ", max_key_len - 1, " bytes")
    }
    if (length(filename_raw) > max_filename_len) {
      stop("Filename '", filenames[i], "' is too long. Maximum length is ", max_filename_len - 1, " bytes")
    }
    
    # Copy data to matrices
    keys_array[1:length(key_raw), i] <- key_raw
    filenames_array[1:length(filename_raw), i] <- filename_raw
  }
  
  ierr <- integer(1)
  
  # Call the generic Fortran subroutine
  result <- .Fortran("create_zip_archive_generic_R",
               zip_raw, as.integer(zip_len),
               as.raw(keys_array), as.integer(max_key_len), as.integer(count),
               as.raw(filenames_array), as.integer(max_filename_len), as.integer(count),
               ierr = ierr)
  
  check_err_code(result$ierr)
  message("Successfully created archive: ", zip_filename)
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
  
  # Flags initialization
  gene_ids_array_valid <- gene_ids_name_valid <- FALSE
  expression_vectors_array_valid <- expression_vectors_name_valid <- FALSE
  gene_to_fam_array_valid <- gene_to_fam_name_valid <- FALSE
  family_ids_array_valid <- family_ids_name_valid <- FALSE
  family_centroids_array_valid <- family_centroids_name_valid <- FALSE
  shift_vectors_array_valid <- shift_vectors_name_valid <- FALSE
  
  ierr <- integer(1)
  if (!is.character(zip_filename)) {
    stop("Type mismatch: zip_filename must be a string")
  }
  
  # Validations
  if (!is.null(gene_ids)) {
    if (length(dim(gene_ids)) == 1 || is.vector(gene_ids)) {
      gene_ids_array_valid <- TRUE
    } else {
      stop(paste("Gene IDs dimensions mismatch: Expected 1 but got", length(dim(gene_ids))))
    }
  }
  if (!is.null(gene_ids_name)) {
    if (is.character(gene_ids_name)) {
      gene_ids_name_valid <- TRUE
    } else {
      stop("Gene IDs name must be a string")
    }
  } else{
    gene_ids_name <- ""
  }
  if (xor(gene_ids_array_valid, gene_ids_name_valid)) {
    if(debug) {message("Gene IDs: Either provide array and filename or none. Skipping.")}
  }
  
  if (!is.null(expression_vectors)) {
    if (length(dim(expression_vectors)) == 2) {
      expression_vectors_array_valid <- TRUE
    } else {
      stop(paste("Expression vectors dim mismatch: Expected 2 but got", length(dim(expression_vectors))))
    }
  }
  if (!is.null(expression_vectors_name)) {
    if (is.character(expression_vectors_name)) {
      expression_vectors_name_valid <- TRUE
    } else {
      stop("Expression vector name must be a string")
    }
  } else {
    expression_vectors_name <- ""
  }
  if (xor(expression_vectors_array_valid, expression_vectors_name_valid)) {
    if(debug) {message("Expression vectors: Either provide array and filename or none. Skipping.")}
  }
  
  if (!is.null(gene_to_fam)) {
    if (length(dim(gene_to_fam)) == 1 || is.vector(gene_to_fam)) {
      gene_to_fam_array_valid <- TRUE
    } else {
      stop(paste("Gene to family dim mismatch: Expected 1 but got", length(dim(gene_to_fam))))
    }
  }
  if (!is.null(gene_to_fam_name)) {
    if (is.character(gene_to_fam_name)) {
      gene_to_fam_name_valid <- TRUE
    } else {
      stop("Gene to family name must be a string")
    }
  } else {
    gene_to_fam_name=""
  }
  if (xor(gene_to_fam_array_valid, gene_to_fam_name_valid)) {
    if(debug) {message("Gene to family: Either provide array and filename or none. Skipping.")}
  }
  
  if (!is.null(family_ids)) {
    if (length(dim(family_ids)) == 1 || is.vector(family_ids)) {
      family_ids_array_valid <- TRUE
    } else {
      stop(paste("Family IDs dim mismatch: Expected 1 but got", length(dim(family_ids))))
    }
  }
  if (!is.null(family_ids_name)) {
    if (is.character(family_ids_name)) {
      family_ids_name_valid <- TRUE
    } else {
      stop("Family IDs name must be a string")
    }
  } else{
    family_ids_name <- ""
  }
  if (xor(family_ids_array_valid, family_ids_name_valid)) {
    if(debug) {message("Family IDs: Either provide array and filename or none. Skipping.")}
  }
  
  if (!is.null(family_centroids)) {
    if (length(dim(family_centroids)) == 2) {
      family_centroids_array_valid <- TRUE
    } else {
      stop(paste("Family centroids dim mismatch: Expected 2 but got", length(dim(family_centroids))))
    }
  }
  if (!is.null(family_centroids_name)) {
    if (is.character(family_centroids_name)) {
      family_centroids_name_valid <- TRUE
    } else {
      stop("Family centroids name must be a string")
    }
  } else{
    family_centroids_name <- ""
  }
  if (xor(family_centroids_array_valid, family_centroids_name_valid)) {
    if(debug) {message("Family centroids: Either provide array and filename or none. Skipping.")}
  }
  
  if (!is.null(shift_vectors)) {
    if (length(dim(shift_vectors)) == 2) {
      shift_vectors_array_valid <- TRUE
    } else {
      stop(paste("Shift vectors dim mismatch: Expected 2 but got", length(dim(shift_vectors))))
    }
  }
  if (!is.null(shift_vectors_name)) {
    if (is.character(shift_vectors_name)) {
      shift_vectors_name_valid <- TRUE
    } else {
      stop("Shift vectors name must be a string")
    }
  } else{
    shift_vectors_name <- ""
  }
  if (xor(shift_vectors_array_valid, shift_vectors_name_valid)) {
    if(debug) {message("Shift vectors: Either provide array and filename or none. Skipping.")}
  }
    
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
  
  if (!is.character(zip_filename) || nchar(zip_filename) == 0) {
    stop("Zip name needs to be a non-empty string")
  }
  ierr <- integer(1)
  
  result <- list(
    gene_ids = NULL,
    expression_vectors = NULL,
    gene_to_fam = NULL,
    family_ids = NULL,
    family_centroids = NULL,
    shift_vectors = NULL
  )

  res <- .Fortran("extract_zip_archive_generic_R", charToRaw(zip_filename), nchar(zip_filename), ierr)

  # Read manifest file
  manifest_path <- "manifest.txt"
  if (!file.exists(manifest_path)) {
    stop("Manifest file not found in archive")
  }
  
  manifest_lines <- readLines(manifest_path)
  file.remove("manifest.txt")
  
  # Parse manifest
  file_mapping <- list()
  for (line in manifest_lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      file_mapping[[parts[1]]] <- parts[2]
    }
  }
  
  # Extract and read files based on parameters
  if (!is.null(gene_ids) && !is.null(file_mapping[["gene_ids"]])) {
    filename <- file_mapping[["gene_ids"]]
    if (file.exists(filename)) {
      result$gene_ids <- tox_deserialize_char_array(filename)
      message(paste("Gene ids extracted from", filename))
      file.remove(filename)
    }
  }
  
  if (!is.null(expression_vectors) && !is.null(file_mapping[["expression"]])) {
    filename <- file_mapping[["expression"]]
    if (file.exists(filename)) {
      result$expression_vectors <- tox_deserialize_real_array(filename)
      message(paste("Expression Vectors extracted from", filename))
      file.remove(filename)
    }
  }
  
  if (!is.null(gene_to_fam) && !is.null(file_mapping[["gene_to_family"]])) {
    filename <- file_mapping[["gene_to_family"]]
    if (file.exists(filename)) {
      result$gene_to_fam <- tox_deserialize_int_array(filename)
      message(paste("Gene to family mapping extracted from", filename))
      file.remove(filename)
    }
  }
  
  if (!is.null(family_ids) && !is.null(file_mapping[["family_ids"]])) {
    filename <- file_mapping[["family_ids"]]
    if (file.exists(filename)) {
      result$family_ids <- tox_deserialize_char_array(filename)
      message(paste("Extracted family IDs from", filename))
      file.remove(filename)
    }
  }
  
  if (!is.null(family_centroids) && !is.null(file_mapping[["family_centroids"]])) {
    filename <- file_mapping[["family_centroids"]]
    if (file.exists(filename)) {
      result$family_centroids <- tox_deserialize_real_array(filename)
      message(paste("Extracted family centroids from", filename))
      file.remove(filename)
    }
  }
  
  if (!is.null(shift_vectors) && !is.null(file_mapping[["shift_vectors"]])) {
    filename <- file_mapping[["shift_vectors"]]
    if (file.exists(filename)) {
      result$shift_vectors <- tox_deserialize_real_array(filename)
      message(paste("Extracted shift vectors from", filename))
      file.remove(filename)
    }
  }
  
  return(result)
}

extract_zip_archive <- function(zip_filename) {
  if (!is.character(zip_filename) || nchar(zip_filename) == 0) {
    stop("Zip name needs to be a non-empty string")
  }
  
  # Fortran Aufruf zur Extraktion
  ierr <- integer(1)
  res <- .Fortran("extract_zip_archive_generic_R", 
                  charToRaw(zip_filename), 
                  nchar(zip_filename), 
                  ierr = ierr)
  
  check_err_code(res$ierr)
  
  # Manifest Datei lesen und Mapping erstellen
  manifest_path <- "manifest.txt"
  if (!file.exists(manifest_path)) {
    stop("Manifest file not found in archive")
  }
  
  manifest_lines <- readLines(manifest_path)
  file_mapping <- list()
  for (line in manifest_lines) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      file_mapping[[parts[1]]] <- parts[2]
    }
  }
  
  return(file_mapping)
}