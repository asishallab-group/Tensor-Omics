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

# --------------------------------
# Atomic Validators
# --------------------------------


validate_is_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

validate_is_integer <- function(x, name = deparse(substitute(x))) {
  if (!is.integer(x)) stop(sprintf("%s must be an integer vector", name))
  invisible(TRUE)
}
validate_is_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1) stop(sprintf("%s must be a scalar", name))
  invisible(TRUE)
}

validate_is_logical <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x)) stop(sprintf("%s must be a logical vector", name))
  invisible(TRUE)
}

validate_is_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) stop(sprintf("%s must be a character vector", name))
  invisible(TRUE)
}

validate_is_vector <- function(x, name = deparse(substitute(x))) {
  if (!is.vector(x)) stop(sprintf("%s must be a vector", name))
  invisible(TRUE)
}

validate_dimensions <- function(x, dim_expected = NULL, name = deparse(substitute(x)), kind = "object") {
  if (is.null(dim_expected)) return(invisible(TRUE))

  dim_expected <- as.integer(dim_expected)
  actual_dim <- dim(x)

  if (is.null(actual_dim)) stop(sprintf("%s must have dimensions", name))
  if (length(actual_dim) != length(dim_expected)) {
    stop(sprintf("%s must be a %d-D %s", name, length(dim_expected), kind))
  }

  for (i in seq_along(dim_expected)) {
    if (!is.na(dim_expected[i]) && actual_dim[i] != dim_expected[i]) {
      stop(sprintf("%s has invalid dimension %d: expected %d, got %d", name, i, dim_expected[i], actual_dim[i]))
    }
  }

  invisible(TRUE)
}

validate_is_matrix <- function(x, name = deparse(substitute(x)), dim = NULL) {
  if (!is.matrix(x)) stop(sprintf("%s must be a matrix", name))
  validate_dimensions(x, dim_expected = dim, name = name, kind = "matrix")
  invisible(TRUE)
}

validate_is_array <- function(x, name = deparse(substitute(x)), dim = NULL) {
  if (!is.array(x)) stop(sprintf("%s must be an array", name))
  validate_dimensions(x, dim_expected = dim, name = name, kind = "array")
  invisible(TRUE)
}

validate_nonempty <- function(x, name = deparse(substitute(x))) {
  if (length(x) == 0) stop(sprintf("%s cannot be empty", name))
  invisible(TRUE)
}

# Backward-compatible alias used by existing wrappers/tests.
validate_nonempty_vector <- function(x, name = deparse(substitute(x))) {
  validate_nonempty(x, name)
}

validate_length_equals <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) stop(sprintf("%s must have length %d", name, n))
  invisible(TRUE)
}

validate_same_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) {
  if (length(x) != length(y)) stop(sprintf("%s and %s must have the same length", name_x, name_y))
  invisible(TRUE)
}

# Backward-compatible alias used by existing wrappers/tests.
validate_equal_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) {
  validate_same_length(x, y, name_x, name_y)
}

validate_length_equals_n <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) stop(sprintf("%s length must equal %d", name, n))
  invisible(TRUE)
}

validate_is_in_range <- function(x, min = -Inf, max = Inf, name = deparse(substitute(x))) {
  if (min > max) stop(sprintf("Invalid range: min (%s) must be <= max (%s)", min, max))
  if (any(x < min, na.rm = TRUE) || any(x > max, na.rm = TRUE)) stop(sprintf("%s must contain values between %s and %s", name, min, max))
  invisible(TRUE)
}

validate_no_na <- function(x, name = deparse(substitute(x))) {
  if (any(is.na(x))) stop(sprintf("%s contains NA values", name))
  invisible(TRUE)
}

validate_no_nan <- function(x, name = deparse(substitute(x))) {
  if (any(is.nan(x))) stop(sprintf("%s contains NaN values", name))
  invisible(TRUE)
}

validate_no_inf <- function(x, name = deparse(substitute(x))) {
  if (any(is.infinite(x))) stop(sprintf("%s contains infinite values", name))
  invisible(TRUE)
}
validate_divisible_length <- function(x, d, name = deparse(substitute(x))) {
  if (d <= 0) stop(sprintf("d must be a positive integer for %s", name))
  if (length(x) %% d != 0) stop(sprintf("Length of %s must be divisible by d", name))
  invisible(TRUE)
}
validate_max_dims <- function(max_dims, name = deparse(substitute(max_dims))) {
  if (!is.numeric(max_dims) || length(max_dims) != 1 || max_dims <= 0) stop(sprintf("%s must be a positive integer", name))
  invisible(TRUE)
}

validate_file_exists <- function(filename, name = deparse(substitute(filename))) {
  if (!file.exists(filename)) stop(sprintf("File does not exist: %s", filename))
  
  invisible(TRUE)
}

validate_filename <- function(filename, name = deparse(substitute(filename))) {
  if (!is.character(filename) || length(filename) != 1) stop(sprintf("%s must be a single character string", name))
  
  invisible(TRUE)
}

validate_non_empty_string <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || nchar(x) == 0L) stop(sprintf("Type mismatch: %s must be a non-empty string", name))
  
  invisible(TRUE)
}



# --------------------------------
# Combined Validators
# --------------------------------

validate_numeric_vector <- function(x, name = deparse(substitute(x))) {
  validate_is_vector(x, name)
  validate_is_numeric(x, name)
  invisible(TRUE)
}
validate_numeric_scalar <- function(x, name = deparse(substitute(x))) {
  validate_is_numeric(x, name)
  validate_is_scalar(x, name)
  invisible(TRUE)
}

validate_numeric_matrix <- function(x, name = deparse(substitute(x)), dim = NULL) {
  validate_is_matrix(x, name = name, dim = dim)
  validate_is_numeric(x, name)
  invisible(TRUE)
}

# Backward-compatible alias used by existing wrappers/tests.
validate_matrix <- function(x, name = deparse(substitute(x)), dim = NULL) {
  validate_numeric_matrix(x, name = name, dim = dim)
}

validate_numeric_array <- function(x, name = deparse(substitute(x)), dim = NULL) {
  validate_is_array(x, name = name, dim = dim)
  validate_is_numeric(x, name)
  invisible(TRUE)
}

validate_character_vector <- function(x, name = deparse(substitute(x))) {
  validate_is_character(x, name)
  invisible(TRUE)
}

validate_character_array <- function(arr, name = deparse(substitute(arr))) {
  validate_is_array(arr, name)
  validate_is_character(arr, name)
  invisible(TRUE)
}

validate_logical_vector <- function(x, name = deparse(substitute(x)),expected_length = NULL) {
  validate_is_logical(x, name)

  if (!is.null(expected_length))
    validate_length_equals(x, expected_length, name)

  invisible(TRUE)
}

validate_integer_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL, bounds = NULL) {
  validate_is_integer(x, name)
  if (!is.null(expected_length)) validate_length_equals(x, expected_length, name)
  if (!is.null(bounds)) {
    lo <- if (!is.null(bounds$min)) bounds$min else -Inf
    hi <- if (!is.null(bounds$max)) bounds$max else Inf
    validate_is_in_range(x, min = lo, max = hi, name = name)
  }
  invisible(TRUE)
}

validate_numeric_matrix_values <- function(x, name = deparse(substitute(x))) {
  validate_numeric_matrix(x, name)
  validate_no_na(x, name)
  validate_no_nan(x, name)
  validate_no_inf(x, name)
  invisible(TRUE)
}

validate_array_or_vector <- function(x, name = deparse(substitute(x))) {
  if (!is.array(x) && !is.vector(x)) stop(sprintf("%s must be an array or vector", name))
  invisible(TRUE)
}

validate_logical_or_index_vector <- function(v, expected_length = NULL, name = deparse(substitute(v))) {
  if (!(is.logical(v) || is.numeric(v) || is.integer(v))) {
    stop(sprintf("%s must be logical or numeric", name))
  }
  if (!is.null(expected_length)) validate_length_equals(v, expected_length, name)
  invisible(TRUE)
}

validate_index_bounds <- function(idx, low = 1, high = Inf, name = deparse(substitute(idx))) {
  if (!is.numeric(idx) && !is.integer(idx)) stop(sprintf("%s must be numeric or integer.", name))
  validate_is_in_range(idx, min = low, max = high, name = name)
  invisible(TRUE)
}


validate_matching_rows <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) {
  validate_is_matrix(A, name_A)
  validate_is_matrix(B, name_B)
  if (nrow(A) != nrow(B)) stop(sprintf("Number of rows in %s must match number of rows in %s.", name_A, name_B))
  invisible(TRUE)
}


validate_matching_cols <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) {
  validate_is_matrix(A, name_A)
  validate_is_matrix(B, name_B)
  if (ncol(A) != ncol(B)) stop(sprintf("Number of columns in %s must match number of columns in %s.", name_A, name_B))
  invisible(TRUE)
}


validate_positive_integer_scalar <- function(x, name = deparse(substitute(x))) {

  validate_numeric_scalar(x, name)
  if (x != as.integer(x)) stop(sprintf("%s must be an integer-valued scalar", name))
  if (as.integer(x) <= 0) stop(sprintf("%s must be a positive integer", name))
  invisible(TRUE)
}

validate_positive_numeric_scalar <- function(x, name = deparse(substitute(x))) {
  validate_numeric_scalar(x, name)
  if (x <= 0) stop(sprintf("%s must be positive", name))
  invisible(TRUE)
}

validate_string_scalar <- function(x, name = deparse(substitute(x))) {

  validate_is_character(x, name)

  if (length(x) != 1) stop(sprintf("%s must be a single string", name))
 invisible(TRUE)
}

validate_index_vector_and_position <- function(ix, position, name = deparse(substitute(ix))) {
  # Validate that ix is an array or vector
  validate_array_or_vector(ix, name)
  
  # Check if index vector is empty
  if (length(ix) == 0) stop(sprintf("Empty %s vector", name))

  # Validate position bounds
  if (position < 1 || position > length(ix)) stop(sprintf("position out of bounds for %s", name))
  invisible(TRUE)
}

validate_group_vectors <- function(group_s, group_c, n_cols,
                                   name_s = deparse(substitute(group_s)),
                                   name_c = deparse(substitute(group_c))) {
  validate_logical_or_index_vector(group_s, name = name_s)
  validate_logical_or_index_vector(group_c, name = name_c)
  validate_same_length(group_s, group_c, name_s, name_c)
  validate_nonempty(group_s, name_s)
  validate_positive_integer_scalar(n_cols, "n_cols")

  if (any(group_s != as.integer(group_s), na.rm = TRUE)) {
    stop(sprintf("%s must contain integer-valued starts", name_s))
  }
  if (any(group_c != as.integer(group_c), na.rm = TRUE)) {
    stop(sprintf("%s must contain integer-valued counts", name_c))
  }

  group_s <- as.integer(group_s)
  group_c <- as.integer(group_c)

  validate_index_bounds(group_s, low = 1L, high = as.integer(n_cols), name = name_s)
  if (any(group_c <= 0L, na.rm = TRUE)) {
    stop(sprintf("%s must contain positive counts", name_c))
  }
  if (any(group_s + group_c - 1L > as.integer(n_cols), na.rm = TRUE)) {
    stop("group_s and group_c define ranges outside input columns")
  }

  invisible(TRUE)
}

# -----------------------------
# Higher-level validators for specific wrappers
# -----------------------------

# Validate gene_to_family vector (integer, length = n_genes, values in 1..n_families)
validate_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) stop(sprintf("%s must be numeric or integer.", name))

  validate_length_equals(gene_to_fam, as.integer(n_genes), name)
  validate_no_na(gene_to_fam, name)
  validate_is_in_range(gene_to_fam, min = 1, max = as.integer(n_families), name = name)
  invisible(TRUE)
}

# Validate gene_to_family for outlier workflow and preserve mapped backend-style error text.
validate_gene_to_family_outliers <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(201)
  }
  validate_length_equals(gene_to_fam, as.integer(n_genes), name)
  validate_no_na(gene_to_fam, name)
  if (any(gene_to_fam < 1 | gene_to_fam > as.integer(n_families), na.rm = TRUE)) {
    check_err_code(201)
  }
  invisible(TRUE)
}

# Validate gene_to_family vector for distance-to-centroid wrappers:
# length must match genes, allow 0 and out-of-range positive values,
# but reject negatives to keep backend semantics used in tests.
validate_gene_to_family_centroid <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) stop(sprintf("%s must be numeric or integer", name))
  if (length(gene_to_fam) != as.integer(n_genes)) stop(sprintf("%s must have equal number of genes", name))
  if (any(gene_to_fam < 0, na.rm = TRUE)) stop(sprintf("%s must be between 0 and %d", name, as.integer(n_families)))
  invisible(TRUE)
}

# Validate ortholog_set vector (logical, length = n_genes)
validate_mode <- function(mode, allowed = c('all', 'ortho', 'orthologs')) {
  if (!mode %in% allowed) stop(sprintf("`mode` must be one of: %s.", paste(allowed, collapse = ", ")))
  invisible(NULL)
}

# Validate inputs for group centroid calculation (expression_vectors, gene_to_family, n_families, ortholog_set, mode).
validate_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {

  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_genes <- ncol(expression_vectors)
  validate_integer_vector(gene_to_family, "gene_to_family", expected_length = n_genes)
  validate_logical_vector(ortholog_set, "ortholog_set", expected_length = n_genes)
  if (!is.numeric(n_families) || length(n_families) != 1 || as.integer(n_families) < 1) stop("`n_families` must be a positive integer scalar.")
  validate_index_bounds(gene_to_family, low = 1, high = as.integer(n_families), name = "gene_to_family")
  validate_mode(mode)
  invisible(NULL)
}
# Validate inputs for mean vector calculation (expression_vectors and gene_indices).
validate_mean_vector_inputs <- function(expression_vectors, gene_indices) {
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_genes <- ncol(expression_vectors)
  validate_integer_vector(gene_indices, name = "gene_indices")
  validate_is_in_range(gene_indices, min = 1, max = n_genes, name = "gene_indices")
  invisible(TRUE)

}

# Validate permutation vectors (1..N, unique).
validate_is_permutation <- function(x, n = length(x), name = deparse(substitute(x))) {
  validate_integer_vector(x, name = name, expected_length = n)
  validate_is_in_range(x, min = 1, max = n, name = name)
  if (length(unique(x)) != n) stop(sprintf("%s must be a valid permutation of 1..%d", name, n))
  invisible(TRUE)
}

# Validate matrix shape for k-means clustering (data_points)
validate_matrix_shape_data_points <- function(data_points, n_dims, n_points, name = "data_points") {
  if (nrow(data_points) != n_dims || ncol(data_points) != n_points) {
    stop(sprintf("%s matrix shape mismatch: expected %%d x %%d, got %%d x %%d", name, n_dims, n_points, nrow(data_points), ncol(data_points)))
  }
  invisible(TRUE)
}

# Validate matrix shape for centroids
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
