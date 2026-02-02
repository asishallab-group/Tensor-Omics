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
# Validation helpers for R wrappers (atomic + composed)
# -----------------------------

# ---- Atomic validators ----
validate_is_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

validate_is_integer <- function(x, name = deparse(substitute(x))) {
  if (!is.integer(x)) stop(sprintf("%s must be an integer vector.", name))
  invisible(TRUE)
}

validate_is_logical <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x)) stop(sprintf("%s must be a logical vector.", name))
  invisible(TRUE)
}

validate_is_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) stop(sprintf("%s must be a character vector", name))
  invisible(TRUE)
}

validate_is_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1) stop(sprintf("%s must be a scalar", name))
  invisible(TRUE)
}

validate_is_vector <- function(x, name = deparse(substitute(x))) {
  if (!is.vector(x)) stop(sprintf("%s must be a vector.", name))
  invisible(TRUE)
}

validate_is_matrix <- function(x, dim = NULL, name = deparse(substitute(x))) {
  if (!is.matrix(x)) stop(sprintf("%s must be a matrix.", name))
  if (!is.null(dim)) {
    if (length(dim) != 2) stop("dim must be length 2 for matrix validation.")
    expected <- dim
    actual <- c(nrow(x), ncol(x))
    if ((!is.na(expected[1]) && actual[1] != expected[1]) || (!is.na(expected[2]) && actual[2] != expected[2])) {
      stop(sprintf("%s matrix shape mismatch: expected %d x %d, got %d x %d", name, expected[1], expected[2], actual[1], actual[2]))
    }
  }
  invisible(TRUE)
}

validate_is_array <- function(x, dim = NULL, name = deparse(substitute(x))) {
  if (!is.array(x)) stop(sprintf("%s must be an array.", name))
  if (!is.null(dim)) {
    actual <- dim(x)
    if (length(actual) != length(dim)) {
      stop(sprintf("%s must have %d dimensions.", name, length(dim)))
    }
    for (i in seq_along(dim)) {
      if (!is.na(dim[i]) && actual[i] != dim[i]) {
        stop(sprintf("%s dimension %d mismatch: expected %s, got %s", name, i, dim[i], actual[i]))
      }
    }
  }
  invisible(TRUE)
}

validate_is_nonempty <- function(x, name = deparse(substitute(x))) {
  if (length(x) == 0) stop(sprintf("%s must not be empty (cannot be empty).", name))
  invisible(TRUE)
}

validate_is_finite <- function(x, name = deparse(substitute(x))) {
  vals <- as.vector(x)
  if (any(is.na(vals))) stop(sprintf("%s contains NA values.", name))
  if (any(is.infinite(vals))) stop(sprintf("%s contains infinite values.", name))
  if (any(is.nan(vals))) stop(sprintf("%s contains NaN values.", name))
  invisible(TRUE)
}

validate_is_in_range <- function(x, min = -Inf, max = Inf, name = deparse(substitute(x))) {
  if (!is.numeric(x) && !is.integer(x)) stop(sprintf("%s must be numeric or integer.", name))
  if (any(!is.finite(as.vector(x)))) stop(sprintf("Invalid input: %s contains NA or non-finite values", name))
  if (any(x < min) || any(x > max)) {
    stop(sprintf("Invalid input: %s indices must be between %s and %s.", name, min, ifelse(is.infinite(max), "Inf", as.character(max))))
  }
  invisible(TRUE)
}

validate_is_same_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) {
  if (length(x) != length(y)) stop(sprintf("%s and %s must have the same length", name_x, name_y))
  invisible(TRUE)
}

validate_is_length_equals <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != as.integer(n)) stop(sprintf("Length of %s must equal number of genes.", name))
  invisible(TRUE)
}

validate_is_divisible_length <- function(x, d, name = deparse(substitute(x))) {
  if (d <= 0) stop(sprintf("d must be a positive integer for %s", name))
  if (length(x) %% d != 0) stop(sprintf("Length of %s must be divisible by d", name))
  invisible(TRUE)
}

validate_is_positive_integer_scalar <- function(x, name = deparse(substitute(x))) {
  if (!(is.numeric(x) || is.integer(x))) stop(sprintf("%s must be numeric or integer", name))
  if (length(x) != 1 || is.na(x) || as.integer(x) <= 0) stop(sprintf("%s must be a positive integer scalar", name))
  invisible(TRUE)
}

validate_is_positive_numeric_scalar <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("Invalid input: %s must be numeric", name))
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) stop(sprintf("Invalid input: %s must be a single positive numeric value", name))
  invisible(TRUE)
}

# ---- Composed validators ----
validate_is_numeric_vector <- function(x, name = deparse(substitute(x))) {
  validate_is_numeric(x, name)
  validate_is_vector(x, name)
  invisible(TRUE)
}

validate_is_numeric_matrix <- function(x, dim = NULL, name = deparse(substitute(x))) {
  validate_is_matrix(x, dim, name)
  if (!(is.numeric(x) || is.integer(x))) {
    stop(sprintf("%s must be a numeric matrix (must be numeric).", name))
  }
  invisible(TRUE)
}

validate_is_numeric_array <- function(x, dim = NULL, name = deparse(substitute(x))) {
  if (!(is.atomic(x) && (is.numeric(x) || is.integer(x)))) {
    stop(sprintf("%s must be a numeric (or integer) vector/array.", name))
  }
  validate_is_nonempty(x, name)
  if (!is.null(dim)) {
    xdim <- dim(x)
    if (is.null(xdim)) xdim <- length(x)
    if (length(dim) == 1) {
      if (!is.na(dim[1]) && xdim != dim[1]) stop(sprintf("%s length mismatch: expected %s, got %s", name, dim[1], xdim))
    } else if (is.null(dim(x)) || length(dim(x)) != length(dim)) {
      stop(sprintf("%s must have %d dimensions.", name, length(dim)))
    } else {
      for (i in seq_along(dim)) {
        if (!is.na(dim[i]) && dim(x)[i] != dim[i]) stop(sprintf("%s dimension %d mismatch: expected %s, got %s", name, i, dim[i], dim(x)[i]))
      }
    }
  }
  invisible(TRUE)
}

validate_is_integer_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL, bounds = NULL) {
  validate_is_integer(x, name)
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  if (!is.null(bounds)) {
    if (!is.null(bounds$min) && any(x < bounds$min)) stop(sprintf("%s contains values < %d.", name, bounds$min))
    if (!is.null(bounds$max) && any(x > bounds$max)) stop(sprintf("%s contains values > %d.", name, bounds$max))
  }
  invisible(TRUE)
}

validate_is_logical_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL) {
  validate_is_logical(x, name)
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("%s must be of length %d.", name, expected_length))
  invisible(TRUE)
}

validate_is_character_vector <- function(x, name = deparse(substitute(x))) {
  validate_is_character(x, name)
  validate_is_vector(x, name)
  invisible(TRUE)
}

validate_is_string_scalar <- function(x, name = deparse(substitute(x))) {
  validate_is_character(x, name)
  validate_is_scalar(x, name)
  invisible(TRUE)
}

validate_is_non_empty_string <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || nchar(x) == 0L) {
    stop(sprintf("Type mismatch: %s must be a non-empty string", name))
  }
  invisible(TRUE)
}

validate_is_filename <- function(filename, name = deparse(substitute(filename))) {
  if (!is.character(filename) || length(filename) != 1) {
    stop(sprintf("%s must be a single character string", name))
  }
  invisible(TRUE)
}

validate_is_file_exists <- function(filename, name = deparse(substitute(filename))) {
  if (!file.exists(filename)) {
    stop(sprintf("File does not exist: %s", filename))
  }
  invisible(TRUE)
}

validate_is_max_dims <- function(max_dims, name = deparse(substitute(max_dims))) {
  if (!is.numeric(max_dims) || length(max_dims) != 1 || max_dims <= 0) {
    stop(sprintf("%s must be a positive integer", name))
  }
  invisible(TRUE)
}

validate_is_array_or_vector <- function(arr, name = deparse(substitute(arr))) {
  if (!is.array(arr) && !is.vector(arr)) {
    stop(sprintf("%s must be an array or vector", name))
  }
  invisible(TRUE)
}

validate_is_character_array <- function(arr, name = deparse(substitute(arr))) {
  if (!is.character(arr)) {
    stop(sprintf("%s must be a character array", name))
  }
  invisible(TRUE)
}

validate_is_group_vectors <- function(group_s, group_c, n_columns) {
  if (!is.integer(group_s) && !is.numeric(group_s)) stop("group_s must be integer (or numeric)")
  if (!is.integer(group_c) && !is.numeric(group_c)) stop("group_c must be integer (or numeric)")
  if (length(group_s) != length(group_c)) stop("group_s and group_c must have the same length")
  invisible(NULL)
}

validate_is_logical_or_index_vector <- function(v, expected_length = NULL, name = deparse(substitute(v))) {
  if (!(is.logical(v) || is.numeric(v) || is.integer(v))) stop(sprintf("%s must be logical or numeric", name))
  if (!is.null(expected_length) && length(v) != expected_length) stop(sprintf("%s length must match expected length %d", name, expected_length))
  invisible(TRUE)
}

validate_is_index_vector <- function(x, n_total, name = deparse(substitute(x))) {
  if (missing(n_total) || !is.numeric(n_total) || length(n_total) != 1 || n_total <= 0) {
    stop("Invalid input: n_total must be a positive integer scalar")
  }
  n_total <- as.integer(n_total)

  if (is.null(x)) return(seq_len(n_total))

  if (is.logical(x)) {
    if (length(x) == 1L) {
      if (!isTRUE(x)) return(integer(0))
      return(seq_len(n_total))
    }
    if (length(x) != n_total) stop(sprintf("Invalid input: %s logical vector length must match expected length %s", name, n_total))
    return(which(x))
  }

  if (!(is.numeric(x) || is.integer(x))) stop(sprintf("Invalid input: %s must be NULL, logical, or numeric index vector", name))
  idx <- as.integer(x)
  if (any(is.na(idx))) stop(sprintf("Invalid input: %s contains NA or non-finite values", name))
  if (length(idx) == 0L) return(integer(0))
  if (any(idx < 1L) || any(idx > n_total)) stop(sprintf("Invalid input: %s indices out of bounds: must be between 1 and %s", name, n_total))
  return(idx)
}

validate_is_index_vector_and_position <- function(ix, position, name = deparse(substitute(ix))) {
  validate_is_array_or_vector(ix, name)
  if (length(ix) == 0) stop(sprintf("Empty %s vector", name))
  if (position < 1 || position > length(ix)) stop(sprintf("position out of bounds for %s", name))
  invisible(TRUE)
}

validate_is_numeric_matrix_values <- function(m, name = deparse(substitute(m))) {
  validate_is_numeric_matrix(m, name = name)
  vals <- as.vector(m)
  na_count <- sum(is.na(vals))
  if (na_count > 0) stop(sprintf("%s contains NA values: %d detected", name, na_count))
  inf_count <- sum(is.infinite(vals))
  if (inf_count > 0) stop(sprintf("%s contains infinite values: %d detected", name, inf_count))
  nan_count <- sum(is.nan(vals))
  if (nan_count > 0) stop(sprintf("%s contains NaN values: %d detected", name, nan_count))
  invisible(NULL)
}

validate_is_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) stop(sprintf("%s must be numeric or integer.", name))
  if (length(gene_to_fam) != as.integer(n_genes)) stop(sprintf("Length of %s must equal number of genes (%d).", name, as.integer(n_genes)))
  if (any(is.na(gene_to_fam))) stop(sprintf("%s contains NA values.", name))
  if (any(gene_to_fam < 0)) stop(sprintf("%s must be between 0 and %d.", name, as.integer(n_families)))
  invisible(TRUE)
}

validate_is_gene_to_centroid <- function(x, name = deparse(substitute(x))) {
  validate_is_numeric(x, name)
  if (any(is.na(x))) stop(sprintf("%s must not contain NA values", name))
  if (any(x < 0L)) stop(sprintf("%s must not contain negative indices", name))
  invisible(TRUE)
}

validate_is_mode <- function(mode, allowed = c('all', 'ortho', 'orthologs')) {
  if (!mode %in% allowed) stop(sprintf("`mode` must be one of: %s.", paste(allowed, collapse = ", ")))
  invisible(NULL)
}

validate_is_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  validate_is_numeric_matrix(expression_vectors, name = "expression_vectors")
  n_genes <- ncol(expression_vectors)
  validate_is_integer_vector(as.integer(gene_to_family), name = "gene_to_family", expected_length = n_genes)
  validate_is_logical_vector(as.logical(ortholog_set), name = "ortholog_set", expected_length = n_genes)
  validate_is_positive_integer_scalar(n_families, "n_families")
  validate_is_mode(mode)
  invisible(NULL)
}

validate_is_mean_vector_inputs <- function(expression_vectors, gene_indices) {
  validate_is_numeric_matrix(expression_vectors, name = "expression_vectors")
  n_genes <- ncol(expression_vectors)
  if (!is.integer(gene_indices)) stop("gene_indices must be an integer vector of 1-based column indices.")
  if (any(gene_indices < 1) || any(gene_indices > n_genes)) stop("gene_indices must contain indices between 1 and ncol(expression_vectors).")
  invisible(TRUE)
}

validate_is_matching_rows <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) {
  validate_is_matrix(A, name = name_A)
  validate_is_matrix(B, name = name_B)
  if (nrow(A) != nrow(B)) stop(sprintf("Number of rows in %s must match number of rows in %s.", name_A, name_B))
  invisible(TRUE)
}

validate_is_loess_smooth_2d_inputs <- function(
    x_ref, y_ref, x_query, indices_used,
    kernel_sigma, kernel_cutoff) {

  validate_is_numeric_vector(x_ref, deparse(substitute(x_ref)))

  A <- y_ref; name <- deparse(substitute(y_ref))
  if (is.null(dim(A))) A <- matrix(A, nrow = 1L)
  validate_is_numeric_matrix(A, name = name)
  if (nrow(A) != 1L) stop(sprintf("%s must have exactly one row.", name))
  if (ncol(A) != length(x_ref)) stop(sprintf("Number of columns in %s must match length(x_ref).", name))

  validate_is_numeric_vector(x_query, deparse(substitute(x_query)))

  if (!is.null(indices_used)) {
    A <- indices_used; name <- deparse(substitute(indices_used))
    validate_is_integer_vector(A, name = name)
    if (any(is.na(A))) stop(sprintf("%s must not contain NA values.", name))
    if (any(A < 1L) || any(A > length(x_ref))) stop(sprintf("%s must contain indices between 1 and length(x_ref).", name))
  }

  validate_is_numeric(kernel_sigma, deparse(substitute(kernel_sigma)))
  validate_is_scalar(kernel_sigma, deparse(substitute(kernel_sigma)))
  validate_is_numeric(kernel_cutoff, deparse(substitute(kernel_cutoff)))
  validate_is_scalar(kernel_cutoff, deparse(substitute(kernel_cutoff)))

  invisible(TRUE)
}

validate_is_matrix_shape_data_points <- function(data_points, n_dims, n_points, name = "data_points") {
  if (nrow(data_points) != n_dims || ncol(data_points) != n_points) {
    stop(sprintf("%s matrix shape mismatch: expected %d x %d, got %d x %d", name, n_dims, n_points, nrow(data_points), ncol(data_points)))
  }
  invisible(TRUE)
}

validate_is_matrix_shape_centroids <- function(centroids, n_dims, n_clusters, name = "centroids") {
  if (nrow(centroids) != n_dims || ncol(centroids) != n_clusters) {
    stop(sprintf("%s matrix shape mismatch: expected %d x %d, got %d x %d", name, n_dims, n_clusters, nrow(centroids), ncol(centroids)))
  }
  invisible(TRUE)
}

validate_is_matrix_shape_factor_centroids <- function(centroids, n_factors, n_clusters, name = "centroids") {
  if (nrow(centroids) != n_factors || ncol(centroids) != n_clusters) {
    stop(sprintf("%s matrix shape mismatch: expected %d x %d, got %d x %d", name, n_factors, n_clusters, nrow(centroids), ncol(centroids)))
  }
  invisible(TRUE)
}

# ---- Compatibility wrappers (temporary; remove after call sites are updated) ----
validate_numeric_vector <- function(x, name = deparse(substitute(x))) validate_is_numeric_vector(x, name)
validate_numeric_matrix <- function(x, name = deparse(substitute(x))) validate_is_numeric_matrix(x, name = name)
validate_numeric_array <- function(x, name = deparse(substitute(x))) validate_is_numeric_array(x, name = name)
validate_nonempty <- function(x, name = deparse(substitute(x))) validate_is_nonempty(x, name)
validate_nonempty_vector <- function(x, name = deparse(substitute(x))) validate_is_nonempty(x, name)
validate_same_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) validate_is_same_length(x, y, name_x, name_y)
validate_equal_length <- function(a, b, name_a = deparse(substitute(a)), name_b = deparse(substitute(b))) validate_is_same_length(a, b, name_a, name_b)
validate_divisible_length <- function(x, d, name = deparse(substitute(x))) validate_is_divisible_length(x, d, name)
validate_positive_integer_scalar <- function(x, name = deparse(substitute(x))) validate_is_positive_integer_scalar(x, name)
validate_positive_numeric_scalar <- function(x, name = deparse(substitute(x))) validate_is_positive_numeric_scalar(x, name)
validate_integer_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL, bounds = NULL) validate_is_integer_vector(x, name, expected_length, bounds)
validate_logical_vector <- function(x, name = deparse(substitute(x)), expected_length = NULL) validate_is_logical_vector(x, name, expected_length)
validate_character_vector <- function(cv, name = deparse(substitute(cv))) validate_is_character_vector(cv, name)
validate_string_scalar <- function(s, name = deparse(substitute(s))) validate_is_string_scalar(s, name)
validate_scalar_character <- function(x, name = deparse(substitute(x))) validate_is_string_scalar(x, name)
validate_logical_or_index_vector <- function(v, expected_length = NULL, name = deparse(substitute(v))) validate_is_logical_or_index_vector(v, expected_length, name)
validate_group_vectors <- function(group_s, group_c, n_columns) validate_is_group_vectors(group_s, group_c, n_columns)
validate_numeric_matrix_values <- function(m, name = deparse(substitute(m))) validate_is_numeric_matrix_values(m, name)
validate_index_vector <- function(x, n_total, name = deparse(substitute(x))) validate_is_index_vector(x, n_total, name)
validate_index_vector_and_position <- function(ix, position, name = deparse(substitute(ix))) validate_is_index_vector_and_position(ix, position, name)
validate_index_bounds <- function(idx, low = 1, high = Inf, name = deparse(substitute(idx))) validate_is_in_range(idx, min = low, max = high, name = name)
validate_matching_rows <- function(A, B, name_A = deparse(substitute(A)), name_B = deparse(substitute(B))) validate_is_matching_rows(A, B, name_A, name_B)
validate_length_equals_n <- function(x, n, name = deparse(substitute(x))) validate_is_length_equals(x, n, name)
validate_filename <- function(filename, name = deparse(substitute(filename))) validate_is_filename(filename, name)
validate_max_dims <- function(max_dims, name = deparse(substitute(max_dims))) validate_is_max_dims(max_dims, name)
validate_array_or_vector <- function(arr, name = deparse(substitute(arr))) validate_is_array_or_vector(arr, name)
validate_character_array <- function(arr, name = deparse(substitute(arr))) validate_is_character_array(arr, name)
validate_file_exists <- function(filename, name = deparse(substitute(filename))) validate_is_file_exists(filename, name)
validate_non_empty_string <- function(x, name = deparse(substitute(x))) validate_is_non_empty_string(x, name)
validate_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) validate_is_gene_to_family(gene_to_fam, n_genes, n_families, name)
validate_gene_to_centroid <- function(x, name = deparse(substitute(x))) validate_is_gene_to_centroid(x, name)
validate_mode <- function(mode, allowed = c('all', 'ortho', 'orthologs')) validate_is_mode(mode, allowed)
validate_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') validate_is_group_centroid_inputs(expression_vectors, gene_to_family, n_families, ortholog_set, mode)
validate_mean_vector_inputs <- function(expression_vectors, gene_indices) validate_is_mean_vector_inputs(expression_vectors, gene_indices)
validate_loess_smooth_2d_inputs <- function(x_ref, y_ref, x_query, indices_used, kernel_sigma, kernel_cutoff) validate_is_loess_smooth_2d_inputs(x_ref, y_ref, x_query, indices_used, kernel_sigma, kernel_cutoff)
validate_matrix_shape_data_points <- function(data_points, n_dims, n_points, name = "data_points") validate_is_matrix_shape_data_points(data_points, n_dims, n_points, name)
validate_matrix_shape_centroids <- function(centroids, n_dims, n_clusters, name = "centroids") validate_is_matrix_shape_centroids(centroids, n_dims, n_clusters, name)
validate_matrix_shape_factor_centroids <- function(centroids, n_factors, n_clusters, name = "centroids") validate_is_matrix_shape_factor_centroids(centroids, n_factors, n_clusters, name)
