source("r/tensoromics_functions_tox_data.R")
source("r/tensoromics_functions.R")

# ---- Example: replicate Fortran test logic in R ----

# Define your file lists (replace with your actual file paths)
file <- c("material/kallisto_sex_data_no_na.tsv")

# Parameters
n_genes <- 88327
n_families <- 15512
gene_len <- 32
family_len <- 32
n_samples <- 67
value_cols <- 2:68

# Allocate result matrices
kallisto_expr <- matrix(0, nrow = n_samples, ncol = n_genes)

# Read gene IDs from first file
cat("Reading gene IDs...\n")
res_gene <- read_gene_ids_from_tsv_file(file, n_genes, gene_len, n_header_rows = 1, gene_col = 1)
gene_ids <- res_gene$gene_ids
cat("ierr (gene ids):", res_gene$ierr, "\n")

# Read 6-replicate files
cat("Reading expression file...\n")
res_expr_6 <- read_expression_vectors_tsv(
  file_list = file,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols,
  delimiter = "\t",
  n_samples = 67
)
cat("ierr :", res_expr_6$ierr, "\n")
kallisto_expr[1:67, ] <- res_expr_6$expression_vectors

# Debug: Check dimensions and types before calling Fortran
cat("DEBUG: gene_ids length:", length(gene_ids), "\n")
cat("DEBUG: kallisto_expr dim:", dim(kallisto_expr), "\n")
cat("DEBUG: value_cols:", value_cols, "\n")
cat("DEBUG: n_samples:", 67, "\n")
cat("DEBUG: gene_len:", gene_len, "\n")
cat("DEBUG: file_len:", max(nchar(file)), "\n")

# Read family mapping
cat("Reading family file...\n")
res_family <- read_orthofinder_file("material/Orthogroups.tsv", gene_ids, n_families, family_len)
gene_family_ids <- res_family$family_ids
gene_to_fam <- res_family$gene_to_fam
cat("ierr (family):", res_family$ierr, "\n")
print(sample(gene_family_ids, 10))

# Filter out genes without family assignments
cat("Filtering unassigned genes...\n")
res_filter <- filter_unassigned_genes(gene_ids, kallisto_expr, gene_to_fam)
cat("ierr (filter):", res_filter$ierr, "\n")
gene_ids <- res_filter$gene_ids
kallisto_expr <- res_filter$expression_vectors
gene_to_fam <- res_filter$gene_to_fam
n_genes_kept <- res_filter$n_genes_kept
cat("n_genes_kept:", n_genes_kept, "\n")

# Data validation
cat("Starting data validation...\n")

# Set parameters for validation
n_genes <- n_genes_kept
n_samples <- nrow(kallisto_expr)

# Validate individual components
cat("Validating gene IDs uniqueness...\n")
validate_string_array_uniqueness(gene_ids, n_genes)

cat("Validating family IDs uniqueness...\n")
validate_string_array_uniqueness(gene_family_ids, n_families)

cat("Validating gene-to-family mapping...\n")
validate_gene_to_family_mapping(gene_to_fam, n_genes, n_families)

cat("Validating expression data...\n")
validate_expression_data(kallisto_expr, n_genes, n_samples)

ortholog_set <- rep(TRUE, n_genes)

family_centroids <- tox_group_centroid(kallisto_expr, gene_to_fam, n_families, 'all')

cat(nrow(family_centroids), "x", ncol(family_centroids), "\n")

cat("Validating family centroids...\n")
validate_family_centroids(family_centroids, n_families, n_samples)

res <- tox_compute_shift_vector_field(kallisto_expr, family_centroids, gene_to_fam)
shift_vectors <- matrix(res, nrow=2*n_samples, ncol=n_genes)

cat("Validating shift vectors...\n")
validate_shift_vectors(shift_vectors, kallisto_expr, family_centroids, gene_to_fam, n_samples)

# Comprehensive validation
cat("Performing comprehensive data validation...\n")
validate_all_data(
  n_genes, n_families, n_samples,
  gene_ids, gene_family_ids,
  gene_to_fam, kallisto_expr,
  family_centroids, shift_vectors
)

cat("All validations passed successfully!\n")

# Additional diagnostic output
cat("\n=== DATA SUMMARY ===\n")
cat("Samples:", n_samples, "\n")
cat("Genes:", n_genes, "\n")
cat("Families:", n_families, "\n")
cat("Expression matrix dimensions:", dim(kallisto_expr), "\n")
cat("Gene ID examples:", head(gene_ids), "\n")
cat("Family ID examples:", head(gene_family_ids), "\n")
cat("Gene-to-family mapping examples:", head(gene_to_fam), "\n")

cat("===Archive tests===\n")
save_tox_data(zip_filename="test_archive_1_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v1_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v1_R.bin", 
                      gene_to_fam=gene_to_fam, gene_to_fam_name="gene_to_fam_v1_R.bin",
                      family_ids=gene_family_ids, family_ids_name="family_ids_v1_R.bin", 
                      family_centroids=family_centroids, family_centroids_name="family_centroids_v1_R.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1_R.bin")

save_tox_data(zip_filename="test_archive_2_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v2_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v2_R.bin")

save_tox_data(zip_filename="test_archive_3_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v3_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v3_R.bin", gene_to_fam=gene_to_fam)

save_tox_data(zip_filename="test_archive_4_R.zip", family_centroids=family_centroids, family_centroids_name="family_centroids_v1_R.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1_R.bin")

save_tox_data(zip_filename="test_archive_5_R.zip")

result_1 <- read_tox_data(zip_filename="test_archive_1_R.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_2 <- read_tox_data(zip_filename="test_archive_1_R.zip",
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_3 <- read_tox_data(zip_filename="test_archive_2_R.zip",
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  shift_vectors=TRUE)
tryCatch({
  result_py <- read_tox_data(zip_filename="test_archive_1_py.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from python archive")

  result_f <- read_tox_data(zip_filename="test_archive_1_f.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from fortran archive")
}, error = function(e) {
  cat("Cross-platform test skipped:", e$message, "\n")
})

cat("=== Testing create_zip_archive directly with dummy arrays ===\n")

# Create dummy arrays
cat("Creating non-standard arrays...\n")

array_3d_int <- 1:60
cat("3D int array length:", length(array_3d_int), "\n")

array_1d_float <- c(0.0, 1.5, -2.3)
cat("1D float array:", array_1d_float, "\n")

array_2d_char <- c("short", "medium_length", "very_long_string_here", "a", "bb", "ccc", "test1", "test2", "test3")
cat("2D char array length:", length(array_2d_char), "\n")

array_1d_bool <- sample(c(TRUE, FALSE), 1000, replace = TRUE)
cat("1D bool array length:", length(array_1d_bool), "\n")

array_complex_real <- matrix(rnorm(16), nrow=4, ncol=4)
array_complex_imag <- matrix(rnorm(16), nrow=4, ncol=4)
cat("Complex real part shape: 4x4\n")
cat("Complex imag part shape: 4x4\n")

cat("Serializing arrays to temporary files...\n")

tox_serialize_int_array(array_3d_int, "temp_3d_int.bin")
tox_serialize_real_array(array_1d_float, "temp_1d_float.bin")
tox_serialize_char_array(array_2d_char, "temp_2d_char.bin")

array_1d_bool_int <- as.integer(array_1d_bool)
tox_serialize_int_array(array_1d_bool_int, "temp_1d_bool.bin")

tox_serialize_real_array(array_complex_real, "temp_complex_real.bin")
tox_serialize_real_array(array_complex_imag, "temp_complex_imag.bin")

cat("Test 1: Direct create_zip_archive call with non-standard arrays\n")
keys <- c(
  "custom_3d_int_data",
  "special_float_array", 
  "string_matrix",
  "boolean_mask",
  "complex_real_part",
  "complex_imag_part"
)

filenames <- c(
  "temp_3d_int.bin",
  "temp_1d_float.bin", 
  "temp_2d_char.bin",
  "temp_1d_bool.bin",
  "temp_complex_real.bin",
  "temp_complex_imag.bin"
)

create_zip_archive("test_non_standard_1.zip", keys, filenames)
cat("Successfully created archive with non-standard arrays\n")

cat("Test 2: Mixed standard and non-standard arrays\n")

tox_serialize_char_array(gene_ids[1:5], "temp_standard_gene_ids.bin")
tox_serialize_real_array(kallisto_expr[1:3, 1:5], "temp_standard_expr.bin")

keys_mixed <- c(keys, "standard_gene_ids", "standard_expression")
filenames_mixed <- c(filenames, "temp_standard_gene_ids.bin", "temp_standard_expr.bin")

create_zip_archive("test_mixed_arrays.zip", keys_mixed, filenames_mixed)
cat("Successfully created archive with mixed standard and non-standard arrays\n")

cat("Test 3: Creating archive with custom keys\n")
create_zip_archive("test_custom_keys.zip",
                  keys = c("experiment_config", "results_summary", "raw_data"),
                  filenames = c("temp_2d_char.bin", "temp_complex_real.bin", "temp_3d_int.bin"))

cat("Test 4: Creating archive with single file\n")
create_zip_archive("test_single_file.zip",
                  keys = c("single_data"),
                  filenames = c("temp_1d_float.bin"))

cat("Test 5: Testing error handling with mismatched arrays\n")
tryCatch({
  create_zip_archive("test_error.zip",
                    keys = c("key1", "key2"),
                    filenames = c("file1.bin"))
}, error = function(e) {
  cat("Expected error caught:", e$message, "\n")
})

cat("\n=== Testing reading and comparing non-standard arrays ===\n")

cat("Test 6a: Reading and verifying test_non_standard_1.zip\n")

file_mapping <- extract_zip_archive("test_non_standard_1.zip")

all_correct <- TRUE

# Integer 3D Array
if (!is.null(file_mapping[["custom_3d_int_data"]])) {
  filename <- file_mapping[["custom_3d_int_data"]]
  if (file.exists(filename)) {
    loaded_3d_int <- tox_deserialize_int_array(filename)
    int_3d_match <- all(loaded_3d_int == array_3d_int)
    cat("custom_3d_int_data arrays match:", int_3d_match, "\n")
    if (!int_3d_match) {
      cat("  Original length:", length(array_3d_int), "Loaded length:", length(loaded_3d_int), "\n")
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("custom_3d_int_data not found in archive\n")
  all_correct <- FALSE
}

# Float 1D Array
if (!is.null(file_mapping[["special_float_array"]])) {
  filename <- file_mapping[["special_float_array"]]
  if (file.exists(filename)) {
    loaded_1d_float <- tox_deserialize_real_array(filename)
    float_match <- TRUE
    for (i in 1:length(array_1d_float)) {
      if (is.nan(array_1d_float[i])) {
        if (!is.nan(loaded_1d_float[i])) float_match <- FALSE
      } else if (is.infinite(array_1d_float[i])) {
        if (!is.infinite(loaded_1d_float[i]) || array_1d_float[i] != loaded_1d_float[i]) float_match <- FALSE
      } else {
        if (abs(array_1d_float[i] - loaded_1d_float[i]) > 1e-10) float_match <- FALSE
      }
    }
    cat("special_float_array arrays match:", float_match, "\n")
    if (!float_match) {
      cat("  Original:", array_1d_float, "\n")
      cat("  Loaded:  ", loaded_1d_float, "\n")
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("special_float_array not found in archive\n")
  all_correct <- FALSE
}

# Character 2D Array
if (!is.null(file_mapping[["string_matrix"]])) {
  filename <- file_mapping[["string_matrix"]]
  if (file.exists(filename)) {
    loaded_2d_char <- tox_deserialize_char_array(filename)
    char_match <- all(loaded_2d_char == array_2d_char)
    cat("string_matrix arrays match:", char_match, "\n")
    if (!char_match) {
      cat("  Original:", array_2d_char, "\n")
      cat("  Loaded:  ", loaded_2d_char, "\n")
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("string_matrix not found in archive\n")
  all_correct <- FALSE
}

# Boolean Array
if (!is.null(file_mapping[["boolean_mask"]])) {
  filename <- file_mapping[["boolean_mask"]]
  if (file.exists(filename)) {
    loaded_bool_int <- tox_deserialize_int_array(filename)
    loaded_bool <- as.logical(loaded_bool_int)
    bool_match <- all(loaded_bool == array_1d_bool)
    cat("boolean_mask arrays match:", bool_match, "\n")
    if (!bool_match) {
      cat("  Original length:", length(array_1d_bool), "Loaded length:", length(loaded_bool), "\n")
      # Zeige nur die ersten 10 Unterschiede
      differences <- which(array_1d_bool != loaded_bool)
      if (length(differences) > 0) {
        cat("  First 10 differences at positions:", head(differences, 10), "\n")
      }
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("boolean_mask not found in archive\n")
  all_correct <- FALSE
}

# Complex Arrays
if (!is.null(file_mapping[["complex_real_part"]])) {
  filename <- file_mapping[["complex_real_part"]]
  if (file.exists(filename)) {
    loaded_complex_real <- tox_deserialize_real_array(filename)

    loaded_complex_real_matrix <- matrix(loaded_complex_real, nrow=4, ncol=4)
    complex_real_match <- all.equal(loaded_complex_real_matrix, array_complex_real, tolerance = 1e-10)
    cat("complex_real_part arrays match:", isTRUE(complex_real_match), "\n")
    if (!isTRUE(complex_real_match)) {
      cat("  Difference:", complex_real_match, "\n")
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("complex_real_part not found in archive\n")
  all_correct <- FALSE
}

if (!is.null(file_mapping[["complex_imag_part"]])) {
  filename <- file_mapping[["complex_imag_part"]]
  if (file.exists(filename)) {
    loaded_complex_imag <- tox_deserialize_real_array(filename)
    loaded_complex_imag_matrix <- matrix(loaded_complex_imag, nrow=4, ncol=4)
    complex_imag_match <- all.equal(loaded_complex_imag_matrix, array_complex_imag, tolerance = 1e-10)
    cat("complex_imag_part arrays match:", isTRUE(complex_imag_match), "\n")
    if (!isTRUE(complex_imag_match)) {
      cat("  Difference:", complex_imag_match, "\n")
      all_correct <- FALSE
    }
    file.remove(filename)
  }
} else {
  cat("complex_imag_part not found in archive\n")
  all_correct <- FALSE
}

cat("Test 6a result: All arrays match =", all_correct, "\n")

cat("Test 6b: Reading and verifying test_mixed_arrays.zip\n")

file_mapping <- extract_zip_archive("test_mixed_arrays.zip")
 
if (!is.null(file_mapping[["standard_gene_ids"]])) {
  filename <- file_mapping[["standard_gene_ids"]]
  if (file.exists(filename)) {
    loaded_genes <- tox_deserialize_char_array(filename)
    genes_match <- all(loaded_genes == gene_ids[1:5])
    cat("standard_gene_ids arrays match:", genes_match, "\n")
    file.remove(filename)
  }
}

if (!is.null(file_mapping[["standard_expression"]])) {
  filename <- file_mapping[["standard_expression"]]
  if (file.exists(filename)) {
    loaded_expr <- tox_deserialize_real_array(filename)
    loaded_expr_matrix <- matrix(loaded_expr, nrow=3, ncol=5)
    expr_match <- all.equal(loaded_expr_matrix, kallisto_expr[1:3, 1:5], tolerance = 1e-10)
    cat("standard_expression arrays match:", isTRUE(expr_match), "\n")
    file.remove(filename)
  }
}

file.remove("manifest.txt")

cat("Checking created archives...\n")
archives <- c("test_non_standard_1.zip", "test_mixed_arrays.zip", 
              "test_custom_keys.zip", "test_single_file.zip")

for (archive in archives) {
  if (file.exists(archive)) {
    cat("Archive created:", archive, "\n")
  } else {
    cat("Archive missing:", archive, "\n")
  }
}

cat("Cleaning up temporary files...\n")
temp_files <- c("temp_3d_int.bin", "temp_1d_float.bin", "temp_2d_char.bin",
                "temp_1d_bool.bin", "temp_complex_real.bin", "temp_complex_imag.bin",
                "temp_standard_gene_ids.bin", "temp_standard_expr.bin")

for (temp_file in temp_files) {
  if (file.exists(temp_file)) {
    file.remove(temp_file)
    cat("Removed:", temp_file, "\n")
  }
}

cat("=== All non-standard array tests completed! ===\n")

cat("=== Reading standard archives ===\n")
result_1 <- read_tox_data(zip_filename="test_archive_1_R.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_2 <- read_tox_data(zip_filename="test_archive_1_R.zip",
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_3 <- read_tox_data(zip_filename="test_archive_2_R.zip",
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  shift_vectors=TRUE)

# Cross-platform Tests
tryCatch({
  result_py <- read_tox_data(zip_filename="test_archive_1_py.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from python archive\n")

  result_f <- read_tox_data(zip_filename="test_archive_1_f.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from fortran archive\n")
}, error = function(e) {
  cat("Cross-platform test skipped:", e$message, "\n")
})

cat("=== All tests completed! ===\n")

cat("=== Archive content validation ===\n")

# List all created archives
archive_files <- list.files(pattern = "\\.zip$")
cat("Created archives:\n")
for (arch in archive_files) {
  file_info <- file.info(arch)
  cat(sprintf("  %s (%.2f KB)\n", arch, file_info$size/1024))
}