#> tox_helper: throw error in error case
check_err_code <- function(ierr) {
  if (ierr == 0) return(invisible(NULL))
  msg <- switch(as.character(ierr),
    # I/O errors
    '101' = "Could not open file.",
    '102' = "Could not read magic number.",
    '103' = "Could not read type code.",
    '104' = "Could not read number of dimensions.",
    '105' = "Could not read array dimensions",
    '106' = "Could not read character length.",
    '107' = "Could not read array data.",
    '112' = "Could not write magic number",
    '113' = "Could not write type code",
    '114' = "Could not write number of dimensions",
    '115' = "Could not write dimensions",
    '116' = "Could not write character length",
    '117' = "Could not write array data",
    # ADD MORE HERE
    
    # FORMAT ERRORS
    '200' = "Invalid format detected.",
    '201' = "Invalid input provided.",
    '202' = "Empty input arrays provided.",
    '203' = "Dimension mismatch detected.",
    '204' = "NaN or Inf found in input data.",
    '205' = "Unsupported data type encountered.",
    '206' = "Array size mismatch detected",

    # MEMORY ERRORS
    '301' = "Memory allocation failed.",
    '302' = "Null pointer reference encountered.",

    # FORTRAN RUNTIME ERRORS
    '5002' = "Fortran runtime error: unit not open / not connected.",

    # Internal errors
    '9001' = "Internal error: unexpected state.",
    '9999' = "Unknown error.",
    paste("Unmapped error code:", ierr)
  )
  stop(msg)
}

#>skip snippets

# -----------------------------
# Validation helpers for R wrappers
# -----------------------------

validate_numeric_vector <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

validate_nonempty <- function(x, name = deparse(substitute(x))) {
  if (length(x) == 0) stop(sprintf("%s must not be empty (cannot be empty).", name))
  invisible(TRUE)
}

validate_same_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) {
  if (length(x) != length(y)) stop(sprintf("%s and %s must have the same length", name_x, name_y))
  invisible(TRUE)
}

validate_divisible_length <- function(x, d, name = deparse(substitute(x))) {
  if (d <= 0) stop(sprintf("d must be a positive integer for %s", name))
  if (length(x) %% d != 0) stop(sprintf("Length of %s must be divisible by d", name))
  invisible(TRUE)
}

validate_matrix <- function(m, name = deparse(substitute(m))) {
  validate_numeric_matrix(m)
}

validate_logical_or_index_vector <- function(v, expected_length = NULL, name = deparse(substitute(v))) {
  if (!(is.logical(v) || is.numeric(v) || is.integer(v))) stop(sprintf("%s must be logical or numeric", name))
  if (!is.null(expected_length) && length(v) != expected_length) stop(sprintf("%s length must match expected length %d", name, expected_length))
  invisible(TRUE)
}

validate_group_vectors <- function(group_s, group_c, n_columns) {
  if (!is.integer(group_s) && !is.numeric(group_s)) stop("group_s must be integer (or numeric)")
  if (!is.integer(group_c) && !is.numeric(group_c)) stop("group_c must be integer (or numeric)")
  if (length(group_s) != length(group_c)) stop("group_s and group_c must have the same length")
  invisible(NULL)
}


# Validate numeric matrix values for NA/Inf/NaN and preserve descriptive messages
validate_numeric_matrix_values <- function(m, name = deparse(substitute(m))) {
  if (!is.matrix(m) && !is.numeric(m)) stop(sprintf("%s must be a numeric matrix", name))
  # coerce but keep original for messages
  vals <- as.vector(m)
  na_count <- sum(is.na(vals))
  if (na_count > 0) stop(sprintf("%s contains NA values: %d detected", name, na_count))
  inf_count <- sum(is.infinite(vals))
  if (inf_count > 0) stop(sprintf("%s contains infinite values: %d detected", name, inf_count))
  nan_count <- sum(is.nan(vals))
  if (nan_count > 0) stop(sprintf("%s contains NaN values: %d detected", name, nan_count))
  invisible(NULL)
}

# Validate single string scalar (e.g. control pattern)
validate_string_scalar <- function(s, name = deparse(substitute(s))) {
  if (!is.character(s) || length(s) != 1) stop(sprintf("%s must be a single string", name))
  invisible(NULL)
}

# Backwards-compatible alias expected by some tests/wrappers
validate_scalar_character <- function(x, name = deparse(substitute(x))) {
  # Reuse validate_string_scalar behaviour
  validate_string_scalar(x)
}

# Validate an index vector used for subsetting operations.
# - `x` may be NULL (meaning all indices), a numeric/integer vector, or logical vector.
# - `n_total` is the length of the full object (must be positive integer).
# - `name` is the parameter name for error messages.
validate_index_vector <- function(x, n_total, name = deparse(substitute(x))) {
  if (missing(n_total) || !is.numeric(n_total) || length(n_total) != 1 || n_total <= 0) {
    stop("Invalid input: n_total must be a positive integer scalar")
  }
  n_total <- as.integer(n_total)

  # NULL => all indices
  if (is.null(x)) return(seq_len(n_total))

  # logical vector handling
  if (is.logical(x)) {
    if (length(x) == 1L) {
      if (!isTRUE(x)) return(integer(0))
      return(seq_len(n_total))
    }
    if (length(x) != n_total) stop(sprintf("Invalid input: %s logical vector length must match expected length %s", name, n_total))
    return(which(x))
  }

  # numeric / integer indexing
  if (!(is.numeric(x) || is.integer(x))) stop(sprintf("Invalid input: %s must be NULL, logical, or numeric index vector", name))
  idx <- as.integer(x)
  if (any(is.na(idx))) stop(sprintf("Invalid input: %s contains NA or non-finite values", name))
  if (length(idx) == 0L) return(integer(0))
  if (any(idx < 1L) || any(idx > n_total)) stop(sprintf("Invalid input: %s indices out of bounds: must be between 1 and %s", name, n_total))
  return(idx)
}

# Validate character vector
validate_character_vector <- function(cv, name = deparse(substitute(cv))) {
  if (!is.character(cv)) stop(sprintf("%s must be a character vector", name))
  invisible(NULL)
 }

# Validate a positive integer scalar (e.g., d)
validate_positive_integer_scalar <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x) && !is.integer(x)) stop(sprintf("%s must be numeric/integer", name))
  if (length(x) != 1) stop(sprintf("%s must be a scalar", name))

  invisible(NULL)
}

# Validate gene_to_family vector (length and index bounds)
validate_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.integer(gene_to_fam) && !is.numeric(gene_to_fam)) stop(sprintf("%s must be integer (or numeric)", name))
  if (length(gene_to_fam) != n_genes) stop(sprintf("Length of %s must equal number of genes (%d)", name, n_genes))
  if (any(gene_to_fam < 0)) stop(sprintf("%s must be between 0 and %d (0 = no family assignment)", name, n_families))
  # Allow indices > n_families here; Fortran / C code handles invalid indices and returns -1 for invalid genes.
  invisible(NULL)
}

validate_gene_to_centroid <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x))stop(sprintf("%s must be numeric", name))
  if (any(is.na(x)))stop(sprintf("%s must not contain NA values", name))
  if (any(x < 0L))stop(sprintf("%s must not contain negative indices", name))

  invisible(TRUE)
 }


validate_positive_integer_scalar <- function(x, name = deparse(substitute(x))) {
  if (!(is.numeric(x) || is.integer(x))) stop(sprintf("%s must be numeric or integer", name))
  if (length(x) != 1 || is.na(x) || as.integer(x) <= 0) stop(sprintf("%s must be a positive integer scalar", name))
  invisible(TRUE)
}

# Validate a positive numeric scalar (non-integer allowed)
validate_positive_numeric_scalar <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("Invalid input: %s must be numeric", name))
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) stop(sprintf("Invalid input: %s must be a single positive numeric value", name))
  invisible(TRUE)
}
validate_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  # Accept integer matrices as numeric too (R often uses integer storage)
  if (!is.matrix(x) || !(is.numeric(x) || is.integer(x))) {
    stop(sprintf("%s must be a numeric matrix (must be numeric).", name))
  }
  invisible(TRUE)
}


# Validate numeric array (vector or array of arbitrary dimensions)
validate_numeric_array <- function(x, name = deparse(substitute(x))) {
  # Accept integer or numeric storage for vectors/arrays of any dimension
  if (!(is.atomic(x) && (is.numeric(x) || is.integer(x)))) {
    stop(sprintf("%s must be a numeric (or integer) vector/array.", name))
  }
  if (length(x) == 0) stop(sprintf("%s must not be empty", name))
  invisible(TRUE)
}

validate_integer_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL, bounds = NULL) {
  if (!is.integer(x)) stop(sprintf("%s must be an integer vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  if (!is.null(bounds)) {
    if (!is.null(bounds$min) && any(x < bounds$min)) stop(sprintf("%s contains values < %d.", name, bounds$min))
    if (!is.null(bounds$max) && any(x > bounds$max)) stop(sprintf("%s contains values > %d.", name, bounds$max))
  }
  invisible(TRUE)
}

validate_logical_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL) {
  if (!is.logical(x)) stop(sprintf("%s must be a logical vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  invisible(TRUE)
}

validate_mode <- function(mode, allowed = c('all', 'ortho')) {
  if (!mode %in% allowed) stop(sprintf("`mode` must be one of: %s.", paste(allowed, collapse = ", ")))
  invisible(NULL)
}

# Higher-level validators for specific wrappers
validate_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  validate_numeric_matrix(expression_vectors)
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)
  validate_integer_vector(as.integer(gene_to_family), expected_length = n_genes)
  validate_logical_vector(as.logical(ortholog_set), expected_length = n_genes)
  if (!is.numeric(n_families) || length(n_families) != 1 || as.integer(n_families) < 1) stop("`n_families` must be a positive integer scalar.")
  validate_mode(mode)
  invisible(NULL)
}

validate_mean_vector_inputs <- function(expression_vectors, gene_indices) {
  validate_numeric_matrix(expression_vectors)
  n_genes <- ncol(expression_vectors)
  if (!is.integer(gene_indices)) stop("gene_indices must be an integer vector of 1-based column indices.")
  if (any(gene_indices < 1) || any(gene_indices > n_genes)) stop("gene_indices must contain indices between 1 and ncol(expression_vectors).")
  invisible(TRUE)
}



validate_equal_length <- function(a, b, name_a = deparse(substitute(a)), name_b = deparse(substitute(b))) {
  if (length(a) != length(b)) stop(sprintf("%s length must match length of %s (same length).", name_a, name_b))
  invisible(TRUE)
}

# Ensure vector is non-empty
validate_nonempty_vector <- function(x, name = deparse(substitute(x))) {
  if (length(x) == 0) stop(sprintf("%s must not be empty (cannot be empty).", name))
  invisible(TRUE)
 
 }


# Validate that two matrices have the same number of rows
validate_matching_rows <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) {
  if (!is.matrix(A) || !is.matrix(B)) stop(sprintf("%s and %s must both be matrices.", name_A, name_B))
  if (nrow(A) != nrow(B)) stop(sprintf("Number of rows in %s must match number of rows in %s.", name_A, name_B))
  invisible(TRUE)
}

# Validate length equals given n (used to match gene_to_fam length to number of genes)
validate_length_equals_n <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != as.integer(n)) stop(sprintf("Length of %s must equal number of genes.", name))
  invisible(TRUE)
}

validate_filename <- function(filename, name = deparse(substitute(filename))) {
  if (!is.character(filename) || length(filename) != 1) {
    stop(sprintf("%s must be a single character string", name))
  }
  invisible(TRUE)
}

#' Validate that max_dims is a positive integer
#'
#' @param max_dims Object to validate
#' @param name Name of the object for error messages
#' @return invisible TRUE if valid, throws error otherwise
validate_max_dims <- function(max_dims, name = deparse(substitute(max_dims))) {
  if (!is.numeric(max_dims) || length(max_dims) != 1 || max_dims <= 0) {
    stop(sprintf("%s must be a positive integer", name))
  }
  invisible(TRUE)
}

#' Validate that arr is an array or vector
#'
#' @param arr Object to validate
#' @param name Name of the object for error messages
#' @return invisible TRUE if valid, throws error otherwise
validate_array_or_vector <- function(arr, name = deparse(substitute(arr))) {
  if (!is.array(arr) && !is.vector(arr)) {
    stop(sprintf("%s must be an array or vector", name))
  }
  invisible(TRUE)
}

#' Validate that arr is a character array
#'
#' @param arr Object to validate
#' @param name Name of the object for error messages
#' @return invisible TRUE if valid, throws error otherwise
validate_character_array <- function(arr, name = deparse(substitute(arr))) {
  if (!is.character(arr)) {
    stop(sprintf("%s must be a character array", name))
  }
  invisible(TRUE)
}

#' Validate that file exists
#'
#' @param filename Path to file to check
#' @param name Name of the object for error messages
#' @return invisible TRUE if valid, throws error otherwise
validate_file_exists <- function(filename, name = deparse(substitute(filename))) {
  if (!file.exists(filename)) {
    stop(sprintf("File does not exist: %s", filename))
  }
  invisible(TRUE)
}
# Ensure string is a non-empty scalar
validate_non_empty_string <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || nchar(x) == 0L) {
    stop(sprintf("Type mismatch: %s must be a non-empty string", name))
  }
  invisible(TRUE)
}


validate_loess_smooth_2d_inputs <- function(
    x_ref, y_ref, x_query, indices_used,
    kernel_sigma, kernel_cutoff) {

  # ----- x_ref -----
  {
    A <- x_ref; name <- deparse(substitute(x_ref))
    if (!is.numeric(A) || !is.vector(A))
      stop(sprintf("%s must be a numeric vector.", name))
  }

  # ----- y_ref -----
  {
    A <- y_ref; name <- deparse(substitute(y_ref))

    # allow input as vector -> convert to 1-row matrix
    if (is.null(dim(A)))
      A <- matrix(A, nrow = 1L)

    if (!is.matrix(A) || !is.numeric(A))
      stop(sprintf("%s must be a numeric matrix.", name))

    if (nrow(A) != 1L)
      stop(sprintf("%s must have exactly one row.", name))

    if (ncol(A) != length(x_ref))
      stop(sprintf("Number of columns in %s must match length(x_ref).", name))
  }

  # ----- x_query -----
  {
    A <- x_query; name <- deparse(substitute(x_query))
    if (!is.numeric(A) || !is.vector(A))
      stop(sprintf("%s must be a numeric vector.", name))
  }

  # ----- indices_used -----
  if (!is.null(indices_used)) {
    A <- indices_used; name <- deparse(substitute(indices_used))

    if (!is.integer(A))
      stop(sprintf("%s must be an integer vector.", name))

    if (any(is.na(A)))
      stop(sprintf("%s must not contain NA values.", name))

    if (any(A < 1L) || any(A > length(x_ref)))
      stop(sprintf("%s must contain indices between 1 and length(x_ref).", name))
  }

  # ----- kernel_sigma -----
  {
    A <- kernel_sigma; name <- deparse(substitute(kernel_sigma))
    if (!is.numeric(A) || length(A) != 1L)
      stop(sprintf("%s must be a single numeric value.", name))
  }

  # ----- kernel_cutoff -----
  {
    A <- kernel_cutoff; name <- deparse(substitute(kernel_cutoff))
    if (!is.numeric(A) || length(A) != 1L)
      stop(sprintf("%s must be a single numeric value.", name))
  }

  invisible(TRUE)
}
validate_index_vector_and_position <- function(ix, position, name = deparse(substitute(ix))) {
  # Validate that ix is an array or vector
  validate_array_or_vector(ix, name)
  
  # Check if index vector is empty
  if (length(ix) == 0) {
    stop(sprintf("Empty %s vector", name))
  }
  
  # Validate position bounds
  if (position < 1 || position > length(ix)) {
    stop(sprintf("position out of bounds for %s", name))
  }
  
  invisible(TRUE)
}
# Validate matrix shape for k-means clustering (data_points)
validate_matrix_shape_data_points <- function(data_points, n_dims, n_points, name = "data_points") {
  if (nrow(data_points) != n_dims || ncol(data_points) != n_points) {
    stop(sprintf("%s matrix shape mismatch: expected %%d x %%d, got %%d x %%d", name, n_dims, n_points, nrow(data_points), ncol(data_points)))
  }
  invisible(TRUE)
}

## Centralized R-side validators (kept inside rcpp/ when r/ must not be changed)
## Messages include specific substrings the tests may assert on.

validate_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("%s must be a numeric matrix (must be numeric).", name))
  }
  invisible(TRUE)
}

validate_integer_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL, bounds = NULL) {
  if (!is.integer(x)) stop(sprintf("%s must be an integer vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  if (!is.null(bounds)) {
    if (!is.null(bounds$min) && any(x < bounds$min)) stop(sprintf("%s contains values < %d.", name, bounds$min))
    if (!is.null(bounds$max) && any(x > bounds$max)) stop(sprintf("%s contains values > %d.", name, bounds$max))
  }
  invisible(TRUE)
}

validate_logical_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL) {
  if (!is.logical(x)) stop(sprintf("%s must be a logical vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  invisible(TRUE)
}
validate_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) stop(sprintf("%s must be numeric or integer.", name))
  if (length(gene_to_fam) != as.integer(n_genes)) stop(sprintf("Length of %s must equal number of genes (%d).", name, as.integer(n_genes)))
  if (any(is.na(gene_to_fam))) stop(sprintf("%s contains NA values.", name))
  if (any(gene_to_fam < 0)) stop(sprintf("%s must be between 0 and %d.", name, as.integer(n_families)))
  invisible(TRUE)
}
validate_mode <- function(mode, allowed = c('all', 'orthologs')) {
  if (!mode %in% allowed) stop(sprintf("`mode` must be one of: %s.", paste(allowed, collapse = ", ")))
  invisible(NULL)
}

# Higher-level validators for specific wrappers
validate_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)
  validate_integer_vector(as.integer(gene_to_family), "gene_to_family", expected_length = n_genes)
  validate_logical_vector(as.logical(ortholog_set), "ortholog_set", expected_length = n_genes)
  if (!is.numeric(n_families) || length(n_families) != 1 || as.integer(n_families) < 1) stop("`n_families` must be a positive integer scalar.")
  validate_mode(mode)
  invisible(NULL)
}

validate_mean_vector_inputs <- function(expression_vectors, gene_indices) {
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_genes <- ncol(expression_vectors)
  if (!is.integer(gene_indices)) stop("gene_indices must be an integer vector of 1-based column indices.")
  if (any(gene_indices < 1) || any(gene_indices > n_genes)) stop("gene_indices must contain indices between 1 and ncol(expression_vectors).")
  invisible(TRUE)
}

validate_equal_length <- function(a, b, name_a = deparse(substitute(a)), name_b = deparse(substitute(b))) {
  if (length(a) != length(b)) stop(sprintf("%s length must match length of %s (same length).", name_a, name_b))
  invisible(TRUE)
}



# Ensure vector is non-empty
validate_nonempty_vector <- function(x, name = deparse(substitute(x))) {
  validate_nonempty(x, name)
}

validate_index_bounds <- function(idx, low = 1, high = Inf, name = deparse(substitute(idx))) {
  if (!is.numeric(idx) && !is.integer(idx)) stop(sprintf("%s must be numeric or integer.", name))
  if (any(idx < low, na.rm = TRUE) || any(idx > high, na.rm = TRUE)) stop(sprintf("Invalid input: %s indices must be between %s and %s.", name, low, ifelse(is.infinite(high), "Inf", as.character(high))))
  invisible(TRUE)
}

# Validate that two matrices have the same number of rows
validate_matching_rows <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) {
  if (!is.matrix(A) || !is.matrix(B)) stop(sprintf("%s and %s must both be matrices.", name_A, name_B))
  if (nrow(A) != nrow(B)) stop(sprintf("Number of rows in %s must match number of rows in %s.", name_A, name_B))
# Validate matrix shape for centroids

  invisible(TRUE)
}

validate_matrix_shape_centroids <- function(centroids, n_dims, n_clusters, name = "centroids") {
  if (nrow(centroids) != n_dims || ncol(centroids) != n_clusters) {
    stop(sprintf("%s matrix shape mismatch: expected %%d x %%d, got %%d x %%d", name, n_dims, n_clusters, nrow(centroids), ncol(centroids)))
  }
  invisible(TRUE)
}

# Validate matrix shape for factor clustering centroids
validate_matrix_shape_factor_centroids <- function(centroids, n_factors, n_clusters, name = "centroids") {
  if (nrow(centroids) != n_factors || ncol(centroids) != n_clusters) {
    stop(sprintf("%s matrix shape mismatch: expected %%d x %%d, got %%d x %%d", name, n_factors, n_clusters, nrow(centroids), ncol(centroids)))
  }
  invisible(TRUE)
}
