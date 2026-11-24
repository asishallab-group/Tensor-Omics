# The `tox_shift_vectors` suite: `compute_shift_vector_field` can be called from R, returns the
# documented type and shape, and raises the documented error, naming the argument. Whether the
# numbers are right is the Fortran suite's job (test/mod_test_tox_shift_vectors/), as the coding
# guide's "Where to Test What" asks; the one value a test checks here is there to show that the
# binding hands the arrays over in the right order and the family ids 1-based.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_TISSUES <- 3L
N_GENES <- 4L
N_FAMILIES <- 2L
# gene_to_fam of a gene in no family
UNASSIGNED <- 0L

# (tissues, genes), every gene different.
.expression_vectors <- function() matrix(seq_len(N_TISSUES * N_GENES), nrow = N_TISSUES, ncol = N_GENES)

# (tissues, families): family f is 10*f in every tissue.
.family_centroids <- function() matrix(rep(c(10, 20), each = N_TISSUES), nrow = N_TISSUES, ncol = N_FAMILIES)

.gene_to_fam <- function() c(1L, 2L, UNASSIGNED, 2L)

# The expected error -- `code` NULL for one the binding raises itself, before the library is
# called -- naming `argument`.
.assert_refused <- function(expr, msg, code = NULL, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(msg, ": nothing was raised", call. = FALSE)
  assert_error(stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              sprintf("%s: expected '%s' named, got '%s'", msg, argument, toString(err$argument)))
}

test_compute_shift_vector_field <- function() {
  result <- compute_shift_vector_field(.expression_vectors(), .family_centroids(), .gene_to_fam())
  assert_true(is.double(result), "expected doubles")
  assert_true(identical(dim(result), c(N_TISSUES, 2L, N_GENES)),
              paste("expected dim", toString(c(N_TISSUES, 2L, N_GENES)), "got", toString(dim(result))))
}

test_compute_shift_vector_field_takes_family_ids_1_based <- function() {
  # the one value: id 2 is the second column of family_centroids, and [, 1, g] the centroid
  result <- compute_shift_vector_field(.expression_vectors(), .family_centroids(), .gene_to_fam())
  assert_true(identical(result[, 1, 2], rep(20, N_TISSUES)),
              paste("centroid of gene 2: got", toString(result[, 1, 2])))
}

test_compute_shift_vector_field_rejects_a_family_id_past_n_families <- function() {
  gene_to_fam <- .gene_to_fam()
  gene_to_fam[N_GENES] <- N_FAMILIES + 1L
  .assert_refused(compute_shift_vector_field(.expression_vectors(), .family_centroids(), gene_to_fam),
                  "a family id past n_families", ERR_INVALID_INPUT, "gene_to_fam")
}

test_compute_shift_vector_field_rejects_a_negative_family_id <- function() {
  gene_to_fam <- .gene_to_fam()
  gene_to_fam[1] <- -1L
  .assert_refused(compute_shift_vector_field(.expression_vectors(), .family_centroids(), gene_to_fam),
                  "a negative family id", ERR_INVALID_INPUT, "gene_to_fam")
}

test_compute_shift_vector_field_rejects_no_genes <- function() {
  # n_genes is read off expression_vectors, so that is the argument named
  .assert_refused(compute_shift_vector_field(matrix(numeric(0), nrow = N_TISSUES, ncol = 0), .family_centroids(),
                                             integer(0)),
                  "no genes", ERR_EMPTY_INPUT, "expression_vectors")
}

test_compute_shift_vector_field_rejects_no_tissues <- function() {
  .assert_refused(compute_shift_vector_field(matrix(numeric(0), nrow = 0, ncol = N_GENES),
                                             matrix(numeric(0), nrow = 0, ncol = N_FAMILIES), .gene_to_fam()),
                  "no tissues", ERR_EMPTY_INPUT, "expression_vectors")
}

test_compute_shift_vector_field_rejects_no_families <- function() {
  # n_families is read off family_centroids
  .assert_refused(compute_shift_vector_field(.expression_vectors(), matrix(numeric(0), nrow = N_TISSUES, ncol = 0),
                                             rep(UNASSIGNED, N_GENES)),
                  "no families", ERR_EMPTY_INPUT, "family_centroids")
}

test_compute_shift_vector_field_rejects_nan_expression <- function() {
  expression_vectors <- .expression_vectors()
  expression_vectors[N_TISSUES, N_GENES] <- NaN
  .assert_refused(compute_shift_vector_field(expression_vectors, .family_centroids(), .gene_to_fam()),
                  "a NaN expression", ERR_NAN_INF, "expression_vectors")
}

test_compute_shift_vector_field_rejects_infinite_centroid <- function() {
  family_centroids <- .family_centroids()
  family_centroids[N_TISSUES, N_FAMILIES] <- Inf
  .assert_refused(compute_shift_vector_field(.expression_vectors(), family_centroids, .gene_to_fam()),
                  "an infinite centroid", ERR_NAN_INF, "family_centroids")
}

test_compute_shift_vector_field_rejects_centroids_of_other_tissues <- function() {
  # the binding's own shape check, before the library is called
  .assert_refused(compute_shift_vector_field(.expression_vectors(), matrix(1, nrow = N_TISSUES + 1L, ncol = N_FAMILIES),
                                             .gene_to_fam()),
                  "centroids over another number of tissues", NULL, "family_centroids")
}

test_compute_shift_vector_field_rejects_a_mapping_of_other_genes <- function() {
  # the binding's own shape check, before the library is called
  .assert_refused(compute_shift_vector_field(.expression_vectors(), .family_centroids(), rep(1L, N_GENES + 1L)),
                  "a mapping of another number of genes", NULL, "gene_to_fam")
}

run_all_tests()
