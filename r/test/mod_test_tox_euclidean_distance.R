# The `tox_euclidean_distance` suite: `euclidean_distance` and `distance_to_centroid` can be called
# from R, return the documented types and shapes, and raise the documented errors, naming the
# argument. Whether the numbers are right is the Fortran suite's job
# (test/mod_test_tox_euclidean_distance/), as the coding guide's "Where to Test What" asks; the one
# value each procedure is checked for here is there to show that the binding hands over the whole
# vectors, and the family ids 1-based.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_TISSUES <- 2L
N_GENES <- 4L
N_FAMILIES <- 2L
# gene_to_fam of a gene in no family
UNASSIGNED <- 0L
# the distance such a gene gets
NO_DISTANCE <- -1

# (tissues, genes): family 1's centroid plus [3, 4], family 2's plus [6, 8], unassigned, family 2's
# centroid itself.
.genes <- function() matrix(c(3, 4, 16, 28, 5, 5, 10, 20), nrow = N_TISSUES, ncol = N_GENES)

# (tissues, families): [0, 0] and [10, 20].
.centroids <- function() matrix(c(0, 0, 10, 20), nrow = N_TISSUES, ncol = N_FAMILIES)

.gene_to_fam <- function() c(1L, 2L, UNASSIGNED, 2L)

# assert_error, and that the error names `argument`. On main, assert_error cannot check the
# argument yet; once it can (its `argument =`), this is assert_error(expr, msg, code, argument = ).
.assert_error_naming <- function(expr, msg, code, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  assert_error(if (!is.null(err)) stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              paste0(msg, ": expected '", argument, "' named, but got ",
                     if (is.null(err$argument)) "none" else sprintf("'%s'", toString(err$argument))))
}

# =====================
# euclidean_distance
# =====================

test_euclidean_distance <- function() {
  result <- euclidean_distance(c(1, 2, 3), c(4, 5, 6))
  assert_true(is.double(result) && length(result) == 1L, "expected a single double")
}

test_euclidean_distance_takes_the_whole_vectors <- function() {
  # the one value: [3, 4] to the origin is 5 only if both elements arrive
  result <- euclidean_distance(c(3, 4), c(0, 0))
  assert_true(identical(result, 5), paste("expected 5, got", result))
}

test_euclidean_distance_rejects_vectors_of_different_lengths <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(euclidean_distance(c(1, 2), c(1, 2, 3)), "vectors of different lengths", NULL, "vec2")
}

test_euclidean_distance_rejects_empty_vectors <- function() {
  # n_elements is read off vec1, so that is the argument named
  .assert_error_naming(euclidean_distance(numeric(0), numeric(0)), "empty vectors", ERR_EMPTY_INPUT, "vec1")
}

test_euclidean_distance_rejects_non_numeric_input <- function() {
  .assert_error_naming(euclidean_distance(c("a", "b"), c(1, 2)), "strings", NULL, "vec1")
}

test_euclidean_distance_rejects_nan <- function() {
  .assert_error_naming(euclidean_distance(c(1, NaN), c(1, 2)), "a NaN in vec1", ERR_NAN_INF, "vec1")
}

test_euclidean_distance_rejects_inf <- function() {
  .assert_error_naming(euclidean_distance(c(1, 2), c(1, -Inf)), "an infinite value in vec2", ERR_NAN_INF, "vec2")
}

# =====================
# distance_to_centroid
# =====================

test_distance_to_centroid <- function() {
  result <- distance_to_centroid(.genes(), .centroids(), .gene_to_fam())
  assert_true(is.double(result), "expected doubles")
  assert_true(is.null(dim(result)) && length(result) == N_GENES,
              paste("expected a vector of", N_GENES, "got dim", toString(dim(result)), "length", length(result)))
}

test_distance_to_centroid_takes_family_ids_1_based <- function() {
  # the one value: id 2 is the second column of centroids, and 0 has no centroid
  result <- distance_to_centroid(.genes(), .centroids(), .gene_to_fam())
  assert_true(identical(result, c(5, 10, NO_DISTANCE, 0)), paste("got", toString(result)))
}

test_distance_to_centroid_rejects_a_family_id_past_n_families <- function() {
  gene_to_fam <- .gene_to_fam()
  gene_to_fam[N_GENES] <- N_FAMILIES + 1L
  .assert_error_naming(distance_to_centroid(.genes(), .centroids(), gene_to_fam),
                       "a family id past n_families", ERR_INVALID_INPUT, "gene_to_fam")
}

test_distance_to_centroid_rejects_a_negative_family_id <- function() {
  gene_to_fam <- .gene_to_fam()
  gene_to_fam[1] <- -1L
  .assert_error_naming(distance_to_centroid(.genes(), .centroids(), gene_to_fam),
                       "a negative family id", ERR_INVALID_INPUT, "gene_to_fam")
}

test_distance_to_centroid_rejects_no_genes <- function() {
  # n_genes is read off genes, so that is the argument named
  .assert_error_naming(distance_to_centroid(matrix(numeric(0), nrow = N_TISSUES, ncol = 0), .centroids(), integer(0)),
                       "no genes", ERR_EMPTY_INPUT, "genes")
}

test_distance_to_centroid_rejects_no_tissues <- function() {
  .assert_error_naming(distance_to_centroid(matrix(numeric(0), nrow = 0, ncol = N_GENES),
                                            matrix(numeric(0), nrow = 0, ncol = N_FAMILIES), .gene_to_fam()),
                       "no tissues", ERR_EMPTY_INPUT, "genes")
}

test_distance_to_centroid_rejects_no_families <- function() {
  # n_families is read off centroids
  .assert_error_naming(distance_to_centroid(.genes(), matrix(numeric(0), nrow = N_TISSUES, ncol = 0),
                                            rep(UNASSIGNED, N_GENES)),
                       "no families", ERR_EMPTY_INPUT, "centroids")
}

test_distance_to_centroid_rejects_nan_expression <- function() {
  genes <- .genes()
  genes[N_TISSUES, N_GENES] <- NaN
  .assert_error_naming(distance_to_centroid(genes, .centroids(), .gene_to_fam()),
                       "a NaN expression", ERR_NAN_INF, "genes")
}

test_distance_to_centroid_rejects_infinite_centroid <- function() {
  centroids <- .centroids()
  centroids[N_TISSUES, N_FAMILIES] <- Inf
  .assert_error_naming(distance_to_centroid(.genes(), centroids, .gene_to_fam()),
                       "an infinite centroid", ERR_NAN_INF, "centroids")
}

test_distance_to_centroid_rejects_centroids_of_other_tissues <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(distance_to_centroid(.genes(), matrix(1, nrow = N_TISSUES + 1L, ncol = N_FAMILIES),
                                            .gene_to_fam()),
                       "centroids over another number of tissues", NULL, "centroids")
}

test_distance_to_centroid_rejects_a_mapping_of_other_genes <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(distance_to_centroid(.genes(), .centroids(), rep(1L, N_GENES + 1L)),
                       "a mapping of another number of genes", NULL, "gene_to_fam")
}

run_all_tests()
