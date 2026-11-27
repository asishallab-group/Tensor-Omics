# The `tox_data_tools` suite: the four readers and filters this module exports
# (`read_gene_ids_from_tsv_file`, `read_expression_vectors_tsv`, `read_orthofinder_file`,
# `get_unassigned_mask`) can be called from R, return the documented types and shapes, and raise
# the documented errors; the binding's own refusals name the argument. Whether the values are
# right is the Fortran suite's job (test/mod_test_tox_data_tools/), as the coding guide's "Where
# to Test What" asks; the one value each reader is checked for here is there to show that the
# binding hands over the column, header and value-column arguments as they were meant.
#
# Every file is written into a temporary directory that the case removes again.
#
# Not called here: `read_gene_ids_from_tsv_file` and `read_expression_vectors_tsv` on a missing
# file. Both close a Fortran unit they never opened, which can close the process's stderr or abort
# it (see the batch's report).

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# the gene_to_fam of a gene in no family
UNASSIGNED <- 0L

# assert_error, and that the error names `argument`. On main, assert_error cannot check the
# argument yet; once it can (its `argument =`), this is assert_error(expr, msg, code, argument = ).
.assert_error_naming <- function(expr, msg, code, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  assert_error(if (!is.null(err)) stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              paste0(msg, ": expected '", argument, "' named, but got ",
                     if (is.null(err$argument)) "none" else sprintf("'%s'", toString(err$argument))))
}

# A new temporary directory; the case removes it with on.exit(unlink(directory, recursive = TRUE))
.scratch <- function() {
  directory <- tempfile("tox_data_tools_")
  dir.create(directory)
  directory
}

# Writes `content` to `name` in `directory`, byte for byte, and returns the path
.write <- function(directory, name, content) {
  path <- file.path(directory, name)
  writeChar(content, path, eos = NULL, useBytes = TRUE)
  path
}

# Three genes, a header, and the values of g2 and g1 swapped against `gene_ids` order
.expression_file <- function(directory) {
  .write(directory, "expression.tsv", "id\ta\tb\ng2\t1.5\t2.5\ng1\t3.5\t4.5\ng3\t5.5\t6.5\n")
}

# Two families: OG1 lists g3 and g1, OG2 lists g2; g4 is in none
.orthogroups_file <- function(directory) {
  .write(directory, "Orthogroups.tsv", "Orthogroup\thuman\nOG1\tg3, g1\nOG2\tg2\n")
}

# =====================
# read_gene_ids_from_tsv_file
# =====================

test_read_gene_ids_from_tsv_file <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tsv <- .write(directory, "genes.tsv", "h\tid\nx\tg1\ny\tg22\n")
  result <- read_gene_ids_from_tsv_file(tsv, 8L, 2L, 1L, 2L)
  assert_true(is.character(result) && is.null(dim(result)), "expected a character vector")
  # the one value: column 2 after one header line, so gene_col and n_header_rows arrive
  assert_true(identical(result, c("g1", "g22")), paste("got", toString(result)))
}

test_read_gene_ids_from_tsv_file_reads_a_file_of_no_gene <- function() {
  # moved from mod_test_tox_data.R's test_an_array_of_no_strings (#221 B3)
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tsv <- .write(directory, "no_genes.tsv", "id\tv1\n")
  ids <- read_gene_ids_from_tsv_file(tsv, 8L, 0L, n_header_rows = 1L, gene_col = 1L)
  assert_true(is.character(ids) && length(ids) == 0, "a file of no genes reads as none")
}

test_read_gene_ids_from_tsv_file_refuses_more_rows_than_ids <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tsv <- .write(directory, "genes.tsv", "id\ng1\ng2\ng3\n")
  assert_error(read_gene_ids_from_tsv_file(tsv, 8L, 2L, 1L, 1L), "three rows for two ids", ERR_SIZE_MISMATCH)
}

test_read_gene_ids_from_tsv_file_refuses_a_negative_length_or_count_by_name <- function() {
  # both size the strings the binding allocates before the file is even opened (#215)
  .assert_error_naming(read_gene_ids_from_tsv_file("absent.test.tsv", -1L, 3L, 1L, 1L),
                       "gene_ids_strlen = -1", NULL, "gene_ids_strlen")
  .assert_error_naming(read_gene_ids_from_tsv_file("absent.test.tsv", 8L, -1L, 1L, 1L),
                       "n_gene_ids_elements = -1", NULL, "n_gene_ids_elements")
}

test_read_gene_ids_from_tsv_file_refuses_wrong_types_by_name <- function() {
  .assert_error_naming(read_gene_ids_from_tsv_file(42, 8L, 3L, 1L, 1L), "a number for a filename", NULL, "filename")
  .assert_error_naming(read_gene_ids_from_tsv_file(c("a.tsv", "b.tsv"), 8L, 3L, 1L, 1L), "two filenames", NULL,
                       "filename")
  .assert_error_naming(read_gene_ids_from_tsv_file("absent.test.tsv", 8L, 3L, "one", 1L), "a string for n_header_rows",
                       NULL, "n_header_rows")
  .assert_error_naming(read_gene_ids_from_tsv_file("absent.test.tsv", 8L, 3L, 1L, 1.5), "a fraction for gene_col",
                       NULL, "gene_col")
}

# =====================
# read_expression_vectors_tsv
# =====================

test_read_expression_vectors_tsv <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  expression <- matrix(0, nrow = 2, ncol = 3)
  result <- read_expression_vectors_tsv(.expression_file(directory), c("g1", "g2", "g3"), expression,
                                        1L, 1L, c(3L, 2L), 1L)
  assert_true(is.double(result) && identical(dim(result), c(2L, 3L)),
              paste("expected a 2 x 3 numeric matrix, got dim", toString(dim(result))))
  # the one value: g1 is the file's second row, and value_cols c(3, 2) puts column 3 first
  assert_true(identical(result[, 1], c(4.5, 3.5)), paste("got", toString(result[, 1])))
}

test_read_expression_vectors_tsv_reads_no_file <- function() {
  # moved from mod_test_tox_data.R's test_an_array_of_no_strings (#221 B3): no file and no gene
  result <- read_expression_vectors_tsv(character(0), character(0), matrix(1, 1, 0), 1L, 1L, 2L, 1L)
  assert_true(identical(dim(result), c(1L, 0L)), paste("got dim", toString(dim(result))))
}

test_read_expression_vectors_tsv_takes_a_delimiter <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  csv <- .write(directory, "expression.csv", "id,a\ng1,2.5\n")
  result <- read_expression_vectors_tsv(csv, "g1", matrix(0, 1, 1), 1L, 1L, 2L, 1L, delimiter = ",")
  assert_true(identical(result[1, 1], 2.5), paste("got", result[1, 1]))
}

test_read_expression_vectors_tsv_refuses_wrong_types_by_name <- function() {
  .assert_error_naming(read_expression_vectors_tsv(character(0), character(0), c(1, 2), 1L, 1L, 2L, 1L),
                       "a vector for expression_vectors", NULL, "expression_vectors")
  .assert_error_naming(read_expression_vectors_tsv(character(0), character(0), matrix("a", 1, 0), 1L, 1L, 2L, 1L),
                       "strings for expression_vectors", NULL, "expression_vectors")
  .assert_error_naming(read_expression_vectors_tsv(1, character(0), matrix(1, 1, 0), 1L, 1L, 2L, 1L),
                       "a number for file_list", NULL, "file_list")
  .assert_error_naming(read_expression_vectors_tsv(character(0), character(0), matrix(1, 1, 0), 1L, 1L, 2.5, 1L),
                       "a fraction in value_cols", NULL, "value_cols")
  .assert_error_naming(read_expression_vectors_tsv(character(0), character(0), matrix(1, 1, 0), 1L, 1L, 2L, 1L,
                                                   delimiter = c(",", ";")),
                       "two delimiters", NULL, "delimiter")
  .assert_error_naming(read_expression_vectors_tsv(character(0), character(0), matrix(1, 1, 0), 1L, 1L, 2L, 1L,
                                                   delimiter = ", "),
                       "a delimiter of two characters", NULL, "delimiter")
}

test_read_expression_vectors_tsv_refuses_no_value_column <- function() {
  assert_error(read_expression_vectors_tsv(character(0), character(0), matrix(1, 1, 0), 1L, 1L, integer(0), 1L),
               "no value column", ERR_INVALID_INPUT)
}

test_read_expression_vectors_tsv_refuses_a_column_per_gene_too_few <- function() {
  assert_error(read_expression_vectors_tsv(character(0), c("g1", "g2"), matrix(1, 1, 1), 1L, 1L, 2L, 1L),
               "one column for two genes", ERR_SIZE_MISMATCH)
}

test_read_expression_vectors_tsv_refuses_nan <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tsv <- .write(directory, "expression.tsv", "id\ta\ng1\tNaN\n")
  assert_error(read_expression_vectors_tsv(tsv, "g1", matrix(0, 1, 1), 1L, 1L, 2L, 1L),
               "a NaN in the file", ERR_INVALID_INPUT)
}

# =====================
# read_orthofinder_file
# =====================

test_read_orthofinder_file <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- read_orthofinder_file(.orthogroups_file(directory), c("g1", "g2", "g3", "g4"), 8L, 2L, 4L)
  assert_true(is.list(result) && setequal(names(result), c("family_ids", "gene_to_fam")),
              paste("got", toString(names(result))))
  assert_true(is.character(result$family_ids) && length(result$family_ids) == 2L, "expected two family ids")
  assert_true(is.integer(result$gene_to_fam) && is.null(dim(result$gene_to_fam)) && length(result$gene_to_fam) == 4L,
              "expected an integer vector of 4")
  # the one value: family indices are 1-based, and a gene in no family gets the sentinel
  assert_true(identical(result$family_ids, c("OG1", "OG2")), paste("got", toString(result$family_ids)))
  assert_true(identical(result$gene_to_fam, c(1L, 2L, 1L, UNASSIGNED)), paste("got", toString(result$gene_to_fam)))
}

test_read_orthofinder_file_reads_a_file_of_no_family <- function() {
  # moved from mod_test_tox_data.R's test_an_array_of_no_strings (#221 B3)
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tsv <- .write(directory, "no_families.tsv", "Orthogroup\tsp1\n")
  families <- read_orthofinder_file(tsv, character(0), 8L, 0L, 0L)
  assert_true(length(families$family_ids) == 0 && length(families$gene_to_fam) == 0,
              "a family file of no family reads as none")
}

test_read_orthofinder_file_refuses_a_missing_file <- function() {
  directory <- .scratch()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  assert_error(read_orthofinder_file(file.path(directory, "missing.tsv"), "g1", 8L, 1L, 1L),
               "a missing file", ERR_FILE_OPEN)
}

test_read_orthofinder_file_refuses_a_negative_length_or_count_by_name <- function() {
  .assert_error_naming(read_orthofinder_file("absent.test.tsv", "g1", -1L, 1L, 1L), "family_ids_strlen = -1", NULL,
                       "family_ids_strlen")
  .assert_error_naming(read_orthofinder_file("absent.test.tsv", "g1", 8L, -1L, 1L), "n_family_ids_elements = -1", NULL,
                       "n_family_ids_elements")
  .assert_error_naming(read_orthofinder_file("absent.test.tsv", "g1", 8L, 1L, -1L), "n_gene_to_fam_elements = -1",
                       NULL, "n_gene_to_fam_elements")
}

test_read_orthofinder_file_refuses_wrong_types_by_name <- function() {
  .assert_error_naming(read_orthofinder_file("absent.test.tsv", 1:2, 8L, 1L, 1L), "numbers for gene_ids", NULL,
                       "gene_ids")
  .assert_error_naming(read_orthofinder_file(NULL, "g1", 8L, 1L, 1L), "NULL for a filename", NULL, "filename")
}

# =====================
# get_unassigned_mask
# =====================

test_get_unassigned_mask <- function() {
  result <- get_unassigned_mask(c(2L, UNASSIGNED, 1L))
  assert_true(is.list(result) && setequal(names(result), c("mask", "n_genes_kept")),
              paste("got", toString(names(result))))
  assert_true(is.logical(result$mask) && is.null(dim(result$mask)) && length(result$mask) == 3L,
              "expected a logical vector of 3")
  assert_true(is.integer(result$n_genes_kept) && length(result$n_genes_kept) == 1L, "expected a single integer")
  # the one value: the sentinel's gene is the one dropped
  assert_true(identical(result$mask, c(TRUE, FALSE, TRUE)) && result$n_genes_kept == 2L,
              paste("got", toString(result$mask), result$n_genes_kept))
}

test_get_unassigned_mask_of_no_gene <- function() {
  result <- get_unassigned_mask(integer(0))
  assert_true(length(result$mask) == 0L && result$n_genes_kept == 0L, "no gene, none kept")
}

test_get_unassigned_mask_refuses_wrong_types_by_name <- function() {
  .assert_error_naming(get_unassigned_mask(c(1.5, 2)), "a fraction", NULL, "gene_to_fam")
  .assert_error_naming(get_unassigned_mask("OG1"), "a string", NULL, "gene_to_fam")
}

run_all_tests()
