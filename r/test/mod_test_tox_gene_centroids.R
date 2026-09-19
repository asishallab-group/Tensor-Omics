# The `tox_gene_centroids` suite: each published procedure can be called from R, returns the
# documented type and shape, and raises the documented error. Whether the numbers are right is the
# Fortran suite's job (test/mod_test_tox_gene_centroids.F90), as the coding guide's "Where to Test
# What" asks.
#
# The one value per procedure picks single columns, so no arithmetic is involved: it only shows
# that gene indices and family ids cross the binding 1-based and that the axes stay along the rows.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_AXES <- 2L
N_GENES <- 4L

# (axes, genes); every column is distinct.
.expression_vectors <- function() {
  rbind(c(1, 2, 3, 4),
        c(10, 20, 30, 40))
}

.assert_double_matrix <- function(result, dims, name) {
  assert_true(is.matrix(result) && is.double(result), sprintf("%s: expected a double matrix", name))
  assert_true(identical(dim(result), as.integer(dims)),
              sprintf("%s: expected dim %s, got %s", name, toString(dims), toString(dim(result))))
}

# -----------------------------------------------------------------------------------------------
# mean_vector
# -----------------------------------------------------------------------------------------------

test_mean_vector <- function() {
  centroid <- mean_vector(.expression_vectors(), 3L)
  assert_true(is.double(centroid) && is.null(dim(centroid)) && length(centroid) == N_AXES,
              "centroid: expected a double vector, one per axis")
  # the one value: a single 1-based index selects that column, the axes down its rows
  assert_true(identical(centroid, c(3, 30)), paste("centroid: got", toString(centroid)))
}

test_mean_vector_accepts_an_empty_selection <- function() {
  # no genes selected is allowed (n_selected_genes >= 0), not an error
  centroid <- mean_vector(.expression_vectors(), integer(0))
  assert_true(is.double(centroid) && length(centroid) == N_AXES, "centroid: expected a double vector, one per axis")
}

test_mean_vector_rejects_an_index_out_of_range <- function() {
  # indices are 1-based: 0 and n_genes + 1 both lie outside
  assert_error(mean_vector(.expression_vectors(), 0L), "gene index 0", ERR_INVALID_INPUT)
  assert_error(mean_vector(.expression_vectors(), N_GENES + 1L), "gene index n_genes + 1", ERR_INVALID_INPUT)
}

test_mean_vector_rejects_more_indices_than_genes <- function() {
  # n_selected_genes, read off gene_indices, may not exceed n_genes
  assert_error(mean_vector(.expression_vectors(), rep(1L, N_GENES + 1L)), "more indices than genes",
               ERR_INVALID_INPUT)
}

test_mean_vector_rejects_no_axes <- function() {
  assert_error(mean_vector(matrix(numeric(0), nrow = 0, ncol = N_GENES), 1L), "no axes", ERR_EMPTY_INPUT)
}

test_mean_vector_rejects_na <- function() {
  # R's missing value is a NaN to the library
  expression_vectors <- .expression_vectors()
  expression_vectors[2, 3] <- NA
  assert_error(mean_vector(expression_vectors, 1L), "a missing expression", ERR_NAN_INF)
}

test_mean_vector_rejects_a_vector <- function() {
  # the binding's own type check, before the library is called
  assert_error(mean_vector(c(1, 2, 3), 1L), "expression_vectors must be a matrix")
}

# -----------------------------------------------------------------------------------------------
# group_centroid_orthologs
# -----------------------------------------------------------------------------------------------

test_group_centroid_orthologs <- function() {
  # gene 2 is unassigned (0); gene 4 is in family 2 but not an ortholog
  gene_to_family <- c(2L, 0L, 1L, 2L)
  ortholog_set <- c(TRUE, TRUE, TRUE, FALSE)
  centroid_matrix <- group_centroid_orthologs(.expression_vectors(), gene_to_family, 2L, ortholog_set)
  .assert_double_matrix(centroid_matrix, c(N_AXES, 2L), "centroid_matrix")
  # the one value: each family keeps exactly one ortholog, so its centroid is that gene's column --
  # family 1 is gene 3, family 2 is gene 1
  assert_true(identical(centroid_matrix, cbind(c(3, 30), c(1, 10))),
              paste("centroid_matrix: got", toString(centroid_matrix)))
}

test_group_centroid_orthologs_rejects_a_family_out_of_range <- function() {
  ortholog_set <- rep(TRUE, N_GENES)
  assert_error(group_centroid_orthologs(.expression_vectors(), c(1L, 3L, 1L, 1L), 2L, ortholog_set),
               "family id n_families + 1", ERR_INVALID_INPUT)
  assert_error(group_centroid_orthologs(.expression_vectors(), c(1L, -1L, 1L, 1L), 2L, ortholog_set),
               "a negative family id", ERR_INVALID_INPUT)
}

test_group_centroid_orthologs_rejects_no_families <- function() {
  assert_error(group_centroid_orthologs(.expression_vectors(), rep(0L, N_GENES), 0L, rep(TRUE, N_GENES)),
               "no families", ERR_EMPTY_INPUT)
}

test_group_centroid_orthologs_rejects_no_genes <- function() {
  assert_error(group_centroid_orthologs(matrix(numeric(0), nrow = N_AXES, ncol = 0), integer(0), 1L, logical(0)),
               "no genes", ERR_EMPTY_INPUT)
}

test_group_centroid_orthologs_rejects_infinity <- function() {
  expression_vectors <- .expression_vectors()
  expression_vectors[1, 4] <- Inf
  assert_error(group_centroid_orthologs(expression_vectors, rep(1L, N_GENES), 1L, rep(TRUE, N_GENES)),
               "an infinite expression", ERR_NAN_INF)
}

test_group_centroid_orthologs_rejects_a_missing_ortholog_flag <- function() {
  # the binding's own NA check, before the library is called
  assert_error(group_centroid_orthologs(.expression_vectors(), rep(1L, N_GENES), 1L, c(TRUE, NA, TRUE, TRUE)),
               "an NA in ortholog_set")
}

test_group_centroid_orthologs_rejects_a_short_ortholog_set <- function() {
  # the binding's own extent check, before the library is called
  assert_error(group_centroid_orthologs(.expression_vectors(), rep(1L, N_GENES), 1L, rep(TRUE, N_GENES - 1L)),
               "ortholog_set shorter than n_genes")
}

# -----------------------------------------------------------------------------------------------
# group_centroid_all
# -----------------------------------------------------------------------------------------------

test_group_centroid_all <- function() {
  # genes 2 and 4 are unassigned (0)
  centroid_matrix <- group_centroid_all(.expression_vectors(), c(2L, 0L, 1L, 0L), 2L)
  .assert_double_matrix(centroid_matrix, c(N_AXES, 2L), "centroid_matrix")
  # the one value: each family holds one gene, so its centroid is that gene's column --
  # family 1 is gene 3, family 2 is gene 1
  assert_true(identical(centroid_matrix, cbind(c(3, 30), c(1, 10))),
              paste("centroid_matrix: got", toString(centroid_matrix)))
}

test_group_centroid_all_rejects_a_family_out_of_range <- function() {
  assert_error(group_centroid_all(.expression_vectors(), c(1L, 1L, 3L, 1L), 2L),
               "family id n_families + 1", ERR_INVALID_INPUT)
}

test_group_centroid_all_rejects_no_families <- function() {
  assert_error(group_centroid_all(.expression_vectors(), rep(0L, N_GENES), 0L), "no families", ERR_EMPTY_INPUT)
}

test_group_centroid_all_rejects_no_axes <- function() {
  assert_error(group_centroid_all(matrix(numeric(0), nrow = 0, ncol = N_GENES), rep(1L, N_GENES), 1L),
               "no axes", ERR_EMPTY_INPUT)
}

test_group_centroid_all_rejects_nan <- function() {
  expression_vectors <- .expression_vectors()
  expression_vectors[2, 1] <- NaN
  assert_error(group_centroid_all(expression_vectors, rep(1L, N_GENES), 1L), "a NaN expression", ERR_NAN_INF)
}

test_group_centroid_all_rejects_a_long_gene_to_family <- function() {
  # the binding's own extent check, before the library is called
  assert_error(group_centroid_all(.expression_vectors(), rep(1L, N_GENES + 1L), 1L),
               "gene_to_family longer than n_genes")
}

run_all_tests()
