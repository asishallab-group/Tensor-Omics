source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# Contract tests for save_flyer_json: the call, the file it writes, and its documented errors.
# R has no JSON parser without extra packages, so the file is compared with the hand-written
# expected file of the Fortran suite, which holds the same data set.

EXPECTED_FILE <- "test/test_files/flyer_expected.json"

reference_data <- function() {
  # R's byte compiler pools constants it finds identical(), and identical(0, -0) is TRUE, so a
  # literal -0.0 may come out as 0.0 or turn a 0.0 negative. The negative zeros are therefore
  # computed at run time.
  negative_zero <- -numeric(1)
  expression_vectors <- matrix(c(1.0, 0.5, -2.0,
                                 0.0, 0.0, 0.25,
                                 1024.0, -0.125, 3.0,
                                 1.5, 2.5, -3.5,
                                 0.75, 8.0, -16.0), nrow = 3L)
  expression_vectors[2, 2] <- negative_zero
  family_centroids <- matrix(c(0.375, 4.25, -8.125,
                               512.5, 0.1875, 0.5,
                               0.0, 0.0, 0.0625), nrow = 3L)
  family_centroids[1, 3] <- negative_zero
  list(
    expression_vectors = expression_vectors,
    family_centroids = family_centroids,
    gene_to_fam = c(2L, 1L, 2L, 0L, 1L),
    is_outlier = c(FALSE, TRUE, FALSE, FALSE, TRUE),
    axis_labels = c("liver", "brain", "heart"),
    family_ids = c("F1", "F2", "F3"),
    gene_ids = c("g1", "g2", "g\u00e8ne", "g4", "g5"),
    gene_species = c("human", "mouse", "human", "fly", "say \"hi\""),
    gene_types = c("ortholog", "paralog", "ortholog", "ortholog", "a\\b")
  )
}

write_flyer <- function(data, path) {
  do.call(save_flyer_json, c(data, list(filename = path)))
}

# Fail unless `expr` raises a tox_error with this code, naming this argument.
assert_tox_error <- function(expr, code, argument, msg) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(msg, ": nothing was raised", call. = FALSE)
  if (!inherits(err, "tox_error") || !identical(as.integer(err$code), as.integer(code)))
    stop(msg, ": expected code ", code, ", got: ", conditionMessage(err), call. = FALSE)
  if (!identical(err$argument, argument))
    stop(msg, ": expected argument '", argument, "', got '", err$argument, "'", call. = FALSE)
  if (!grepl(sprintf("'%s'", argument), conditionMessage(err), fixed = TRUE))
    stop(msg, ": the message does not name '", argument, "'", call. = FALSE)
  invisible(TRUE)
}

test_reference_file_matches_expected <- function() {
  path <- tempfile(fileext = ".json")
  on.exit(if (file.exists(path)) file.remove(path))
  write_flyer(reference_data(), path)

  written <- readBin(path, "raw", n = file.size(path))
  expected <- readBin(EXPECTED_FILE, "raw", n = file.size(EXPECTED_FILE))
  assert_true(identical(written, expected), "the file differs from the expected file, byte for byte")
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  assert_true(identical(lines, readLines(EXPECTED_FILE, encoding = "UTF-8", warn = FALSE)),
              "the file differs from the expected file, line for line")
  assert_equal_int(length(lines), 1L, "the document is a single line")
}

test_documented_errors_name_the_argument <- function() {
  path <- tempfile(fileext = ".json")
  on.exit(if (file.exists(path)) file.remove(path))

  data <- reference_data(); data$expression_vectors[2, 3] <- NaN
  assert_tox_error(write_flyer(data, path), ERR_NAN_INF, "expression_vectors", "NaN")

  data <- reference_data(); data$family_centroids[1, 2] <- Inf
  assert_tox_error(write_flyer(data, path), ERR_NAN_INF, "family_centroids", "Inf")

  data <- reference_data(); data$gene_to_fam <- c(2L, 1L, 4L, 0L, 1L)
  assert_tox_error(write_flyer(data, path), ERR_INVALID_INPUT, "gene_to_fam",
                   "family index above the number of families")

  data <- reference_data(); data$axis_labels <- c("liver", "brain", "liver")
  assert_tox_error(write_flyer(data, path), ERR_INVALID_INPUT, "axis_labels", "repeated axis label")

  data <- reference_data(); data$gene_ids <- c("g1", "g2", "g1", "g4", "g5")
  assert_tox_error(write_flyer(data, path), ERR_INVALID_INPUT, "gene_ids", "repeated gene id")

  data <- reference_data(); data$family_ids <- c("F1", "", "F3")
  assert_tox_error(write_flyer(data, path), ERR_INVALID_INPUT, "family_ids", "empty family id")

  assert_false(file.exists(path), "a refused call creates no file")
}

test_existing_file_is_refused <- function() {
  path <- tempfile(fileext = ".json")
  on.exit(if (file.exists(path)) file.remove(path))
  writeLines("keep me", path)

  assert_tox_error(write_flyer(reference_data(), path), ERR_FILE_OPEN, "filename", "existing file")
  assert_true(identical(readLines(path), "keep me"), "the existing file is untouched")
}

test_shape_mismatch_is_refused_before_the_call <- function() {
  path <- tempfile(fileext = ".json")
  data <- reference_data(); data$gene_ids <- c("g1", "g2", "g3", "g4")
  assert_error(write_flyer(data, path), "gene_ids shorter than the genes")
  assert_false(file.exists(path), "a refused call creates no file")
}

run_all_tests()
