# The pure-R helpers in r/data_prep.R. No Fortran stands behind them, so unlike the binding suites
# this one checks values: R is the only place they can be tested.
#
# Both helpers treat ROWS as genes (`n_genes <- nrow(input_matrix)`, NA handling per row), whatever
# their roxygen says; these cases pin that behaviour down.

source("r/test_helpers.R")
source("r/data_prep.R")

# the cleaner narrates every step on stdout
.quietly <- function(expr) {
  capture.output(value <- expr)
  value
}

# rows are genes: one with a gap, one all zero, one with an infinity, one clean
.dirty <- function() {
  rbind(c(1, NA, 3),
        c(0, 0, 0),
        c(2, Inf, 4),
        c(5, 6, 7))
}

test_diagnose_counts_every_kind_of_problem <- function() {
  diagnostics <- tox_diagnose_data_quality(rbind(c(1, NA, 3), c(0, 0, 0), c(-1, Inf, NaN), c(2, 3, 4)),
                                           show_details = FALSE)
  assert_true(identical(unname(diagnostics$dimensions), c(4L, 3L, 12L)), "dimensions")
  problems <- diagnostics$problems
  # is.na() is TRUE for NaN as well, so the NA count includes it
  assert_true(problems$na_count == 2L, "na_count")
  assert_true(problems$inf_count == 1L, "inf_count")
  assert_true(problems$nan_count == 1L, "nan_count")
  assert_true(problems$zero_count == 3L, "zero_count")
  assert_true(problems$negative_count == 1L, "negative_count")
  genes <- lapply(diagnostics$problematic_genes, unname)
  assert_true(identical(genes$genes_with_na, c(1L, 3L)), "genes_with_na")
  assert_true(identical(genes$genes_with_inf, 3L), "genes_with_inf")
  assert_true(identical(genes$genes_with_nan, 3L), "genes_with_nan")
  assert_true(identical(genes$genes_all_zero, 2L), "genes_all_zero")
  # a matrix with problems gets no summary statistics
  assert_true(is.na(diagnostics$statistics$mean_val), "statistics of a dirty matrix")
}

test_diagnose_summarizes_a_clean_matrix <- function() {
  statistics <- tox_diagnose_data_quality(matrix(c(1, 2, 3, 4), nrow = 2), show_details = FALSE)$statistics
  assert_equal_numeric(c(statistics$min_val, statistics$max_val, statistics$mean_val), c(1, 4, 2.5),
                       msg = "min, max, mean")
}

test_diagnose_prints_its_report_only_when_asked <- function() {
  assert_true(length(capture.output(tox_diagnose_data_quality(.dirty()))) > 0L, "the report is printed by default")
  assert_true(length(capture.output(tox_diagnose_data_quality(.dirty(), show_details = FALSE))) == 0L,
              "show_details = FALSE prints nothing")
}

test_clean_removes_genes_with_gaps <- function() {
  cleaned <- .quietly(tox_clean_data_for_normalization(.dirty(), na_strategy = "remove_genes"))
  # the infinity became a gap, and the all-zero gene went too
  assert_true(identical(cleaned, rbind(c(5, 6, 7))), "remove_genes")
}

test_clean_removes_samples_with_gaps <- function() {
  cleaned <- .quietly(tox_clean_data_for_normalization(.dirty(), na_strategy = "remove_samples"))
  assert_true(identical(cleaned, rbind(c(1, 3), c(2, 4), c(5, 7))), "remove_samples")
}

test_clean_imputes_zero <- function() {
  cleaned <- .quietly(tox_clean_data_for_normalization(.dirty(), na_strategy = "impute_zero"))
  assert_true(identical(cleaned, rbind(c(1, 0, 3), c(2, 0, 4), c(5, 6, 7))), "impute_zero")
}

test_clean_imputes_the_gene_mean <- function() {
  cleaned <- .quietly(tox_clean_data_for_normalization(.dirty(), na_strategy = "impute_mean"))
  assert_true(identical(cleaned, rbind(c(1, 2, 3), c(2, 3, 4), c(5, 6, 7))), "impute_mean")
}

test_clean_smart_impute_drops_mostly_empty_genes <- function() {
  # the first gene is two thirds gaps, past the 50% limit; the second has one gap to impute
  cleaned <- .quietly(tox_clean_data_for_normalization(rbind(c(NA, NA, 1), c(1, NA, 3), c(4, 5, 6)),
                                                       na_strategy = "smart_impute"))
  assert_true(identical(cleaned, rbind(c(1, 2, 3), c(4, 5, 6))), "smart_impute")
}

test_clean_keeps_all_zero_genes_when_asked <- function() {
  cleaned <- .quietly(tox_clean_data_for_normalization(.dirty(), remove_all_zero_genes = FALSE,
                                                       na_strategy = "impute_zero"))
  assert_true(identical(cleaned[2, ], c(0, 0, 0)), "the all-zero gene is kept")
}

test_clean_zeroes_small_values_only_when_asked <- function() {
  input <- rbind(c(1, 2, 3), c(0.5, 2, 3))
  kept <- .quietly(tox_clean_data_for_normalization(input, min_expression_threshold = 1.5))
  assert_true(identical(kept, input), "small values survive by default")
  zeroed <- .quietly(tox_clean_data_for_normalization(input, min_expression_threshold = 1.5,
                                                      convert_small_to_zero = TRUE))
  assert_true(identical(zeroed, rbind(c(0, 2, 3), c(0, 2, 3))), "values below the threshold become zero")
}

run_all_tests()
