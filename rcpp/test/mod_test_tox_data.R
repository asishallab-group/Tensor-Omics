source("rcpp/load_tensor_omics.R")
source("rcpp/test_helpers.R")

# ---- Example: replicate Fortran test logic in R ----

test_1_read_expression_vectors_tsv <- function() {
  # Define your file lists (replace with your actual file paths)
  file <- c("material/kallisto_sex_data_no_na.tsv")

  # Parameters
  n_genes <- 88327
  gene_len <- 32
  n_samples <- 67
  value_cols <- 2:68

  # Allocate result matrices
  kallisto_expr <<- matrix(0, nrow = n_samples, ncol = n_genes)

  # one file, and the string width comes before the count now
  gene_ids <<- read_gene_ids_from_tsv_file(file[1], gene_len, n_genes,
                                           n_header_rows = 1, gene_col = 1)

  # expression_vectors is filled in place and returned
  kallisto_expr <<- read_expression_vectors_tsv(
    file_list = file,
    gene_ids = gene_ids,
    expression_vectors = kallisto_expr,
    n_header_rows = 1,
    gene_col = 1,
    value_cols = as.integer(value_cols),
    start_row = 1L,
    delimiter = "\t"
  )

  # Debug: Check dimensions and types before calling Fortran
  assert_equal_int(length(gene_ids), 88327L, "Read gene id count from tsv doesn't match")
}

test_2_read_orthofinder_file <- function() {
  # Read family mapping
  n_families <<- 15512
  family_len <- 32
  res_family <- read_orthofinder_file("material/Orthogroups.tsv", gene_ids,
                                      family_len, n_families, length(gene_ids))
  gene_family_ids <<- res_family$family_ids
  gene_to_fam <<- res_family$gene_to_fam
  assert_true(identical(
    gene_family_ids[gene_to_fam[seq_len(10)]],
    c("OG0000047", "OG0000047", "OG0000047", "OG0000047", "OG0002145", "OG0000230", "OG0000127", "OG0000127", "OG0000127", "OG0000127")
  ), "Family mapping mismatch")
}

test_3_filter_unassigned_genes <- function() {
  gene_to_fam[1] = 0L # unset assignment for first gene, so at least one is unassigned
  # Filter out genes without family assignments

  res_filter <- get_unassigned_mask(gene_to_fam)

  assert_true(identical(gene_to_fam > 0, res_filter$mask), "Filtered mask mismatch")
  gene_ids <<- gene_ids[res_filter$mask]
  kallisto_expr <<- kallisto_expr[, res_filter$mask]
  gene_to_fam <<- gene_to_fam[res_filter$mask]
  n_genes_kept <<- res_filter$n_genes_kept
}

test_4_validation <- function() {
  # Set parameters for validation
  n_genes <- n_genes_kept
  n_samples <- nrow(kallisto_expr)

  # Validate individual components
  force(validate_string_array_uniqueness(gene_ids))

  force(validate_string_array_uniqueness(gene_family_ids))

  force(validate_gene_to_family_mapping(gene_to_fam, n_families))

  force(validate_expression_data(kallisto_expr, check_non_negative = TRUE))

  ortholog_set <- rep(TRUE, n_genes)
  family_centroids <- group_centroid(kallisto_expr, gene_to_fam, n_families,
                                     mode = 'group_all',
                                     ortholog_set = rep(TRUE, length(gene_to_fam)))

  force(validate_family_centroids(family_centroids))

  # the field is a natural (n_samples, 2, n_genes) array now; the validators and the
  # archive still use the flattened (2 * n_samples, n_genes) layout
  res <- compute_shift_vector_field(kallisto_expr, family_centroids, gene_to_fam)
  shift_vectors <- rbind(res[, 1, ], res[, 2, ])

  force(validate_shift_vectors(shift_vectors, kallisto_expr, family_centroids,
                               gene_to_fam, n_samples))

  force(validate_all_data(
    n_genes, n_families, n_samples,
    gene_ids, gene_family_ids,
    gene_to_fam, kallisto_expr,
    family_centroids, shift_vectors
  ))

  invisible()
}

test_5_zipping <- function() {
  # Small synthetic data rather than the 88k-gene pipeline above: this is about the
  # archive round-tripping, and the earlier tests' locals are not all in scope here.
  n_genes <- 4L; n_samples <- 3L; n_families <- 2L
  ids <- sprintf("gene%02d", seq_len(n_genes))
  fams <- sprintf("OG%04d", seq_len(n_families))
  expr <- matrix(as.double(seq_len(n_samples * n_genes)), nrow = n_samples)
  g2f <- as.integer(c(1, 1, 2, 2))
  centroids <- matrix(as.double(seq_len(n_samples * n_families)), nrow = n_samples)
  shifts <- matrix(as.double(seq_len(2 * n_samples * n_genes)), nrow = 2 * n_samples)

  archive <- "archive_roundtrip_R.test.zip"
  if (file.exists(archive)) file.remove(archive)

  save_tox_data(zip_filename = archive,
                gene_ids = ids, gene_ids_file = "gene_ids_R.test.bin",
                expression = expr, expression_file = "expression_R.test.bin",
                gene_to_family = g2f, gene_to_family_file = "gene_to_family_R.test.bin",
                family_ids = fams, family_ids_file = "family_ids_R.test.bin",
                family_centroids = centroids, family_centroids_file = "centroids_R.test.bin",
                shift_vectors = shifts, shift_vectors_file = "shift_vectors_R.test.bin")
  assert_true(file.exists(archive), "the archive should have been written")

  # get_tox_data_dims reports every member's shape, and read_tox_data_into sizes its
  # outputs from it, so the whole archive comes back in one call
  dims <- get_tox_data_dims(archive)
  assert_true(dims$n_gene_ids == n_genes, "gene id count from the archive header")
  assert_true(dims$n_expression_rows == n_samples, "expression rows from the header")
  assert_true(dims$n_family_ids == n_families, "family id count from the header")

  back <- read_tox_data_into(archive)
  assert_true(all(back$gene_ids == ids), "gene ids should round-trip")
  assert_true(all(back$family_ids == fams), "family ids should round-trip")
  assert_true(all(back$gene_to_family == g2f), "gene-to-family should round-trip")
  assert_equal_numeric(as.vector(back$expression), as.vector(expr), 0.0,
                       "expression should round-trip")
  assert_equal_numeric(as.vector(back$family_centroids), as.vector(centroids), 0.0,
                       "centroids should round-trip")
  assert_equal_numeric(as.vector(back$shift_vectors), as.vector(shifts), 0.0,
                       "shift vectors should round-trip")

  # a partial archive: the members left out come back empty rather than erroring
  partial <- "archive_partial_R.test.zip"
  if (file.exists(partial)) file.remove(partial)
  save_tox_data(zip_filename = partial,
                gene_ids = ids, gene_ids_file = "gene_ids_p_R.test.bin")
  partial_dims <- get_tox_data_dims(partial)
  assert_true(partial_dims$n_gene_ids == n_genes, "the member that is there")
  assert_true(partial_dims$n_expression_cols == 0L, "the member that is not")

  assert_error(read_tox_data_into("no_such_archive_R.test.zip"),
               "Expected an error for an archive that does not exist")

  # create_zip_archive directly, with arrays that are not part of the tox data set
  serialize_int_helper(array(seq_len(60), dim = c(5, 4, 3)), "temp_3d_int_R.test.bin")
  serialize_real_helper(c(0.0, 1.5, -2.3), "temp_1d_float_R.test.bin")
  serialize_char_helper(c("short", "medium_length", "a"), "temp_char_R.test.bin")

  keys <- c("custom_3d_int", "special_floats", "strings")
  files <- c("temp_3d_int_R.test.bin", "temp_1d_float_R.test.bin", "temp_char_R.test.bin")
  bundle <- "non_standard_R.test.zip"
  if (file.exists(bundle)) file.remove(bundle)
  create_zip_archive(bundle, keys, files)
  assert_true(file.exists(bundle), "the bundle should have been written")

  # the manifest maps the keys to the member files: a count, the two string widths, then
  # alternating key and filename
  extracted <- tempfile("tox_manifest_")
  dir.create(extracted)
  utils::unzip(bundle, exdir = extracted)
  manifest <- readLines(file.path(extracted, "manifest.txt"))
  n_entries <- as.integer(manifest[1])
  assert_true(n_entries == length(keys), "the manifest should list every key")
  body <- manifest[-seq_len(3)]
  mapped <- trimws(body[seq(1, 2 * n_entries, by = 2)])
  assert_true(setequal(mapped, keys), "the manifest keys should be the ones written")

  # and the members deserialize back to what went in
  member <- file.path(extracted, "temp_3d_int_R.test.bin")
  meta <- get_array_metadata(member, 8L)
  shape <- meta$dims_out[seq_len(meta$ndims)]
  restored <- deserialize_int_helper(shape, member)
  dim(restored) <- shape
  assert_true(all(dim(restored) == c(5, 4, 3)), "the 3d array keeps its shape")
  assert_true(all(as.vector(restored) == seq_len(60)), "and its values")

  assert_error(create_zip_archive("error_R.test.zip", c("key1", "key2"), c("only_one.bin")),
               "Expected an error for mismatched keys and filenames")

  for (f in c(files, archive, partial, bundle)) if (file.exists(f)) file.remove(f)
}

run_all_tests()
