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

validate_numeric_array <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x) || !is.array(x)) stop(sprintf("%s must be numeric array", name))
  invisible(TRUE)
}

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

validate_numeric_matrix <- function(m, name = deparse(substitute(m))) {
  if (!is.matrix(m)) stop(sprintf("%s must be a matrix", name))
  if (!is.numeric(m)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

# Backwards-compatible alias: some wrappers call `validate_matrix`
validate_matrix <- function(m, name = deparse(substitute(m))) {
  validate_numeric_matrix(m, name)
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

validate_gene_to_family <- function(gene_to_fam, n_genes, n_families,
                                    name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {stop(sprintf("%s must be numeric or integer.", name))}
  if (length(gene_to_fam) != as.integer(n_genes)) {stop(sprintf("Length of %s must equal number of genes (%d).", name, as.integer(n_genes)))}
  invisible(TRUE)
}


## Small helpers for messages used in tests
validate_length_equals_n <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) stop(sprintf("%s length must equal %d", name, n))
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
  invisible(TRUE)
}
