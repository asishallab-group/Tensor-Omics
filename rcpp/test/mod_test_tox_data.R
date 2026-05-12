
source("rcpp/tensoromics_functions_tox_data.R")
source("rcpp/test_helpers.R")


make_dataset <- function() {
  gene_ids <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")
  expression <- matrix(
    c(10.0, 12.0, 14.0, 16.0,
      20.0, 22.0, 24.0, 26.0,
      30.0, 32.0, 34.0, 36.0),
    nrow = 3,
    ncol = 4,
    byrow = FALSE,
    dimnames = list(paste0("S", 1:3), gene_ids)
  )
  gene_to_fam <- c(1L, 1L, 2L, 2L) # Two genes in each family
  list(gene_ids = gene_ids, expression = expression, gene_to_fam = gene_to_fam)
}

test_removes <- function() {
  dataset <- make_dataset()
  # Simulate a zero family entry for gene 2
  gene_to_fam_with_zero <- c(1L, 0L, 2L, 2L)
  filtered <- filter_unassigned_genes(dataset$gene_ids, dataset$expression, gene_to_fam_with_zero)
  assert_true(filtered$n_genes_kept == 3L)
  assert_true(identical(filtered$gene_ids, dataset$gene_ids[c(1, 3, 4)]))
  assert_true(ncol(filtered$expression_vectors) == 3L)
}

test_TSV <- function() {
  dataset <- make_dataset()
  tsv_path <- file.path(tempdir(), "tox_data_small.tsv")
  lines <- c(
    "gene\ts1\ts2\ts3",
    sprintf("%s\t%.1f\t%.1f\t%.1f", dataset$gene_ids[1], dataset$expression[1, 1], dataset$expression[2, 1], dataset$expression[3, 1]),
    sprintf("%s\t%.1f\t%.1f\t%.1f", dataset$gene_ids[2], dataset$expression[1, 2], dataset$expression[2, 2], dataset$expression[3, 2]),
    sprintf("%s\t%.1f\t%.1f\t%.1f", dataset$gene_ids[3], dataset$expression[1, 3], dataset$expression[2, 3], dataset$expression[3, 3]),
    sprintf("%s\t%.1f\t%.1f\t%.1f", dataset$gene_ids[4], dataset$expression[1, 4], dataset$expression[2, 4], dataset$expression[3, 4])
  )
  writeLines(lines, tsv_path)

  gene_ids_res <- read_gene_ids_from_tsv_file(tsv_path, 4L, 16L, 1L, 1L)
  gene_ids <- gene_ids_res$gene_ids
  assert_true(identical(trimws(gene_ids), dataset$gene_ids))

  expr_res <- read_expression_vectors_tsv(
    file_list = tsv_path,
    gene_ids = gene_ids,
    n_samples = 3L,
    n_header_rows = 1L,
    gene_col = 1L,
    value_cols = as.integer(2:4)
  )

  expr <- expr_res$expression_vectors
  # Ignore dimnames for comparison
  assert_true(isTRUE(all.equal(unname(expr), unname(dataset$expression), tolerance = 1e-8)))

  unlink(tsv_path)
}


test_validate_string_array_uniqueness <- function() {
  assert_error(validate_string_array_uniqueness(c("A", "A")), "Expected error for duplicates")
}

test_create_zip_archive <- function() {
  assert_error(create_zip_archive(file.path(tempdir(), "bad.test.zip"), keys = "only", filenames = character()), "Expected error for not existing archive")
}

run_all_tests()
