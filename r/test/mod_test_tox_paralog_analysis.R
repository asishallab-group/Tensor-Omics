# The `tox_paralog_analysis` suite: `detect_neofunctionalization`, `detect_dosage_effect`,
# `detect_subfunctionalization`, the two `filter_paralogs_by_pattern_*` filters,
# `calc_paralog_norms` and `calc_paralog_subsets_bitmasks_size` can be called from R, return the documented types and
# shapes, and raise the documented errors, naming the argument. Whether the subsets and flags are
# right is the Fortran suite's job (test/mod_test_tox_paralog_analysis/), as the coding guide's
# "Where to Test What" asks. The few values checked here are the binding's own: genes are bits in
# index order, from 1, a bit mask word with bit 31 set is NA_integer_ and round-trips as one, and
# the binding sizes the work array and cuts max_subset_size itself.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# gene_to_fam of a gene in no family (detect_neofunctionalization only)
UNASSIGNED <- 0L
# the word whose only set bit is bit 31, gene 32 of a mask: 0x80000000, which R holds as NA
ONLY_BIT_31 <- NA_integer_

# Ancestors (axes, families) [0.5, 0.5] and [-0.5, 0]; three genes, in families 2, 1 and none;
# thresholds 0.25.
.ancestors <- function() matrix(c(0.5, 0.5, -0.5, 0), nrow = 2L)
.neofunc_genes <- function() matrix(c(0.5, 0.5, 0.5, 0, -1, 1), nrow = 2L)
.gene_to_fam <- function() c(2L, 1L, UNASSIGNED)
.thresholds <- function() c(0.25, 0.25)

# Three genes, each 0.75 times the ancestor [1, 1], all kept (the word 7): every pair is 1.5 times
# the ancestor, a dosage effect.
.ancestor <- function() c(1, 1)
.dosage_genes <- function() matrix(0.75, nrow = 2L, ncol = 3L)
ALL_THREE <- 7L

# Ancestor [2, 0], three genes [1, 0], all kept, with their norms and ascending order: every pair
# rebuilds the ancestor.
.subfunc_ancestor <- function() c(2, 0)
.subfunc_genes <- function() matrix(c(1, 0, 1, 0, 1, 0), nrow = 2L)
.norms <- function() c(1, 1, 1)
.ascending <- function() 1:3

# Angles [0, 1, 0], one family.
.angles <- function() c(0, 1, 0)
.one_family <- function() rep(1L, 3L)

# assert_error, and that the error names `argument`. On main, assert_error cannot check the
# argument yet; once it can (its `argument =`, #242), this is assert_error(expr, msg, code, argument = ).
.assert_error_naming <- function(expr, msg, code, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  assert_error(if (!is.null(err)) stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              paste0(msg, ": expected '", argument, "' named, but got ",
                     if (is.null(err$argument)) "none" else sprintf("'%s'", toString(err$argument))))
}

# The words of the results a detection found, in increasing order.
.found_words <- function(result) sort(result$paralog_subsets_bitmasks[1L, seq_len(result$n_results)])

# =====================
# detect_neofunctionalization
# =====================

test_detect_neofunctionalization <- function() {
  result <- detect_neofunctionalization(.ancestors(), .neofunc_genes(), .gene_to_fam(), .thresholds())
  assert_true(is.logical(result), "expected logicals")
  assert_true(identical(dim(result), c(3L, 2L)), paste("expected (genes, axes) = (3, 2), got", toString(dim(result))))
}

test_detect_neofunctionalization_takes_family_ids_1_based <- function() {
  # the one value: id 2 is the second column of ancestors, and 0 has none. Gene 1 is family 1's
  # ancestor exactly, so it is flagged only if compared with the second.
  result <- detect_neofunctionalization(.ancestors(), .neofunc_genes(), .gene_to_fam(), .thresholds())
  expected <- matrix(c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE), nrow = 3L)
  assert_true(identical(result, expected), paste("got", toString(result)))
}

test_detect_neofunctionalization_rejects_a_family_id_past_n_families <- function() {
  gene_to_fam <- .gene_to_fam()
  gene_to_fam[3L] <- 3L
  .assert_error_naming(detect_neofunctionalization(.ancestors(), .neofunc_genes(), gene_to_fam, .thresholds()),
                       "a family id past n_families", ERR_INVALID_INPUT, "gene_to_fam")
}

test_detect_neofunctionalization_rejects_genes_out_of_range <- function() {
  genes <- .neofunc_genes()
  genes[2L, 3L] <- 2
  .assert_error_naming(detect_neofunctionalization(.ancestors(), genes, .gene_to_fam(), .thresholds()),
                       "a gene component above 1", ERR_INVALID_INPUT, "genes")
}

test_detect_neofunctionalization_rejects_nan_ancestor <- function() {
  ancestors <- .ancestors()
  ancestors[2L, 2L] <- NaN
  .assert_error_naming(detect_neofunctionalization(ancestors, .neofunc_genes(), .gene_to_fam(), .thresholds()),
                       "a NaN in ancestors", ERR_NAN_INF, "ancestors")
}

test_detect_neofunctionalization_rejects_infinite_threshold <- function() {
  .assert_error_naming(detect_neofunctionalization(.ancestors(), .neofunc_genes(), .gene_to_fam(), c(0.25, Inf)),
                       "an infinite threshold", ERR_NAN_INF, "thresholds")
}

test_detect_neofunctionalization_rejects_no_families <- function() {
  # n_families is read off ancestors
  .assert_error_naming(detect_neofunctionalization(matrix(numeric(0), nrow = 2L, ncol = 0L), .neofunc_genes(),
                                                   rep(UNASSIGNED, 3L), .thresholds()),
                       "no families", ERR_EMPTY_INPUT, "ancestors")
}

test_detect_neofunctionalization_rejects_no_genes <- function() {
  # n_genes is read off genes
  .assert_error_naming(detect_neofunctionalization(.ancestors(), matrix(numeric(0), nrow = 2L, ncol = 0L), integer(0),
                                                   .thresholds()),
                       "no genes", ERR_EMPTY_INPUT, "genes")
}

test_detect_neofunctionalization_rejects_mismatched_shapes <- function() {
  # the binding's own shape checks, before the library is called
  .assert_error_naming(detect_neofunctionalization(.ancestors(), matrix(0, nrow = 3L, ncol = 3L), .gene_to_fam(),
                                                   .thresholds()),
                       "genes over another number of axes", NULL, "genes")
  .assert_error_naming(detect_neofunctionalization(.ancestors(), .neofunc_genes(), .gene_to_fam()[1:2], .thresholds()),
                       "a mapping of another number of genes", NULL, "gene_to_fam")
  .assert_error_naming(detect_neofunctionalization(.ancestors(), .neofunc_genes(), .gene_to_fam(), 0.25),
                       "thresholds of another number of axes", NULL, "thresholds")
}

# =====================
# detect_dosage_effect
# =====================

test_detect_dosage_effect <- function() {
  result <- detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, 3L)
  assert_true(is.list(result) && identical(sort(names(result)), c("n_results", "paralog_subsets_bitmasks")),
              paste("expected a list of n_results and paralog_subsets_bitmasks, got", toString(names(result))))
  assert_true(is.integer(result$n_results) && length(result$n_results) == 1L, "expected a single integer n_results")
  # the whole work array comes back, as many subsets as the sizing asks for
  work_array_size <- calc_paralog_subsets_bitmasks_size(3L, 3L, ALL_THREE)$work_array_size
  assert_true(is.integer(result$paralog_subsets_bitmasks), "expected integer words")
  assert_true(identical(dim(result$paralog_subsets_bitmasks), c(1L, work_array_size)),
              paste("expected (n_mask_chunks, work_array_size), got", toString(dim(result$paralog_subsets_bitmasks))))
}

test_detect_dosage_effect_numbers_genes_from_bit_0 <- function() {
  # the one value: the three pairs are genes {1, 2}, {1, 3} and {2, 3}, the words 3, 5 and 6
  found <- .found_words(detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, 3L))
  assert_true(identical(found, c(3L, 5L, 6L)), paste("got", toString(found)))
}

test_detect_dosage_effect_cuts_max_subset_size <- function() {
  # the binding sizes the work array itself, and the sizing cuts an over-large max_subset_size
  # to the kept genes, so a search size of 2e9 asks for no more than three genes do
  result <- detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, 2000000000L)
  assert_true(identical(dim(result$paralog_subsets_bitmasks), c(1L, 3L)),
              paste("got", toString(dim(result$paralog_subsets_bitmasks))))
  assert_true(identical(result$n_results, 3L), paste("got", result$n_results))
}

test_a_gene_alone_in_bit_31_survives_the_chain <- function() {
  # Kept from the old suite. Gene 32 is bit 31 of word 1. Kept alone in that word, the word is
  # 0x80000000, which R holds as NA_integer_ (issue #206): every step of filter -> calc ->
  # detect takes it back.
  n_genes <- 33L
  angles <- rep(1.0, n_genes)
  angles[c(32L, 33L)] <- 0.1
  gene_to_fam <- rep(1L, n_genes)

  family_bitmasks <- filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1L, gene_to_fam, 2L)
  kept <- family_bitmasks[, 1L]
  assert_true(identical(kept, c(ONLY_BIT_31, 1L)), "genes 32 and 33: the NA word, then bit 0 of word 2")

  sizing <- calc_paralog_subsets_bitmasks_size(2L, n_genes, kept)
  assert_true(identical(sizing$max_subset_size, 2L), "two genes are kept, so a subset holds two at most")

  # only the two kept genes together exceed the ancestor's magnitude by 10%
  genes <- matrix(0.6, nrow = 2L, ncol = n_genes)
  found <- detect_dosage_effect(c(1.0, 1.0), genes, kept, sizing$max_subset_size, gain_gamma = 0.1)
  assert_true(identical(found$n_results, 1L), "the pair {32, 33} is the one dosage subset")
  assert_true(identical(found$paralog_subsets_bitmasks[, 1L], c(ONLY_BIT_31, 1L)),
              "the subset comes back as the same words, the NA word included")

  # an NA in an ordinary integer argument is still refused: it is no family index
  .assert_error_naming(filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1L, c(NA_integer_, gene_to_fam[-1L]), 2L),
                       "gene_to_fam = NA", NULL, "gene_to_fam")
  assert_true(inherits(tryCatch(filter_paralogs_by_pattern_dosage_effect(
                angles, 0.5, 1L, c(NA_integer_, gene_to_fam[-1L]), 2L), error = function(e) e), "tox_na_error"),
              "gene_to_fam = NA is refused as NA")
}

test_detect_dosage_effect_with_fewer_than_two_kept_genes <- function() {
  # Regression: a family whose filter kept no gene, or one, holds no pair, and the detection has
  # nothing to find: no result, n_results 0. The binding sizes the work array at 0 subsets and
  # passes that on, which the library used to refuse as ERR_EMPTY_INPUT naming n_paralog_subsets
  # -- an argument the caller cannot pass.
  for (kept in c(0L, 2L)) {
    result <- detect_dosage_effect(.ancestor(), .dosage_genes(), kept, 3L)
    assert_true(identical(result$n_results, 0L), paste("the word", kept, ": got", result$n_results))
  }
}

test_detect_dosage_effect_rejects_a_negative_max_subset_size <- function() {
  .assert_error_naming(detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, -1L),
                       "max_subset_size -1", ERR_INVALID_INPUT, "max_subset_size")
}

test_detect_dosage_effect_rejects_a_mask_too_short_for_the_genes <- function() {
  # 40 genes need two words; the sizing runs first and blames the genes it was handed
  .assert_error_naming(detect_dosage_effect(.ancestor(), matrix(1, nrow = 2L, ncol = 40L), ALL_THREE, 3L),
                       "one word for 40 genes", ERR_INVALID_INPUT, "genes")
}

test_detect_dosage_effect_rejects_a_gain_of_zero <- function() {
  .assert_error_naming(detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, 3L, gain_gamma = 0),
                       "gain_gamma 0", ERR_INVALID_INPUT, "gain_gamma")
}

test_detect_dosage_effect_rejects_an_angle_above_pi <- function() {
  .assert_error_naming(detect_dosage_effect(.ancestor(), .dosage_genes(), ALL_THREE, 3L, max_angle = 4),
                       "max_angle 4", ERR_INVALID_INPUT, "max_angle")
}

test_detect_dosage_effect_rejects_nan_ancestor <- function() {
  .assert_error_naming(detect_dosage_effect(c(1, NaN), .dosage_genes(), ALL_THREE, 3L),
                       "a NaN in ancestor", ERR_NAN_INF, "ancestor")
}

test_detect_dosage_effect_rejects_infinite_genes <- function() {
  genes <- .dosage_genes()
  genes[2L, 3L] <- -Inf
  .assert_error_naming(detect_dosage_effect(.ancestor(), genes, ALL_THREE, 3L),
                       "an infinite gene", ERR_NAN_INF, "genes")
}

test_detect_dosage_effect_rejects_a_zero_ancestor <- function() {
  .assert_error_naming(detect_dosage_effect(c(0, 0), .dosage_genes(), ALL_THREE, 3L),
                       "an ancestor of zeros", ERR_INVALID_INPUT, "ancestor")
}

test_detect_dosage_effect_rejects_no_dimensions <- function() {
  # n_dims is read off ancestor
  .assert_error_naming(detect_dosage_effect(numeric(0), matrix(numeric(0), nrow = 0L, ncol = 3L), ALL_THREE, 3L),
                       "no dimensions", ERR_EMPTY_INPUT, "ancestor")
}

test_detect_dosage_effect_rejects_genes_of_other_dimensions <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(detect_dosage_effect(c(1, 1, 1), .dosage_genes(), ALL_THREE, 3L),
                       "genes over another number of dimensions", NULL, "genes")
}

# =====================
# detect_subfunctionalization
# =====================

test_detect_subfunctionalization <- function() {
  result <- detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, .norms(), .ascending())
  assert_true(is.list(result) && identical(sort(names(result)), c("n_results", "paralog_subsets_bitmasks")),
              paste("expected a list of n_results and paralog_subsets_bitmasks, got", toString(names(result))))
  assert_true(is.integer(result$n_results) && length(result$n_results) == 1L, "expected a single integer n_results")
  work_array_size <- calc_paralog_subsets_bitmasks_size(3L, 3L, ALL_THREE)$work_array_size
  assert_true(is.integer(result$paralog_subsets_bitmasks), "expected integer words")
  assert_true(identical(dim(result$paralog_subsets_bitmasks), c(1L, work_array_size)),
              paste("expected (n_mask_chunks, work_array_size), got", toString(dim(result$paralog_subsets_bitmasks))))
}

test_detect_subfunctionalization_numbers_genes_from_bit_0 <- function() {
  # the one value: any pair of the three [1, 0] genes rebuilds [2, 0], the words 3, 5 and 6
  found <- .found_words(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, .norms(),
                                                    .ascending()))
  assert_true(identical(found, c(3L, 5L, 6L)), paste("got", toString(found)))
}

test_detect_subfunctionalization_with_fewer_than_two_kept_genes <- function() {
  # Regression: as test_detect_dosage_effect_with_fewer_than_two_kept_genes, which used to be
  # refused as ERR_EMPTY_INPUT naming n_paralog_subsets. No result, n_results 0.
  for (kept in c(0L, 2L)) {
    result <- detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), kept, 3L, 0, .norms(), .ascending())
    assert_true(identical(result$n_results, 0L), paste("the word", kept, ": got", result$n_results))
  }
}

test_detect_subfunctionalization_rejects_a_negative_threshold <- function() {
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, -1, .norms(),
                                                   .ascending()),
                       "rdi_threshold -1", ERR_INVALID_INPUT, "rdi_threshold")
}

test_detect_subfunctionalization_rejects_a_negative_norm <- function() {
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, c(1, 1, -1),
                                                   .ascending()),
                       "a norm of -1", ERR_INVALID_INPUT, "paralog_norms")
}

test_detect_subfunctionalization_rejects_a_permutation_entry_of_0 <- function() {
  # the entries are 1-based gene indices
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, .norms(),
                                                   c(0L, 2L, 3L)),
                       "a permutation entry 0", ERR_INVALID_INPUT, "sorted_paralog_norms_perm")
}

test_detect_subfunctionalization_rejects_a_nan_norm <- function() {
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, c(1, 1, NaN),
                                                   .ascending()),
                       "a NaN norm", ERR_NAN_INF, "paralog_norms")
}

test_detect_subfunctionalization_rejects_mismatched_shapes <- function() {
  # the binding's own shape checks, before the library is called
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, c(1, 1),
                                                   .ascending()),
                       "norms of another number of genes", NULL, "paralog_norms")
  .assert_error_naming(detect_subfunctionalization(.subfunc_ancestor(), .subfunc_genes(), ALL_THREE, 3L, 0, .norms(), 1:2),
                       "a permutation of another number of genes", NULL, "sorted_paralog_norms_perm")
  .assert_error_naming(detect_subfunctionalization(c(2, 0, 0), .subfunc_genes(), ALL_THREE, 3L, 0, .norms(), .ascending()),
                       "genes over another number of dimensions", NULL, "genes")
}

# =====================
# calc_paralog_norms
# =====================

test_calc_paralog_norms <- function() {
  result <- calc_paralog_norms(matrix(c(3, 4, 0, 0.5), nrow = 2L))
  assert_true(is.double(result), "expected doubles")
  assert_true(is.null(dim(result)) && length(result) == 2L,
              paste("expected a vector of 2, got dim", toString(dim(result)), "length", length(result)))
}

test_calc_paralog_norms_takes_genes_as_columns <- function() {
  # the one value: the columns [3, 4] and [0, 0.5] are 5 and 0.5 long
  result <- calc_paralog_norms(matrix(c(3, 4, 0, 0.5), nrow = 2L))
  assert_true(identical(result, c(5, 0.5)), paste("got", toString(result)))
}

test_calc_paralog_norms_rejects_nan_genes <- function() {
  .assert_error_naming(calc_paralog_norms(matrix(c(1, 1, 1, NaN), nrow = 2L)), "a NaN gene", ERR_NAN_INF, "genes")
}

test_calc_paralog_norms_rejects_no_genes <- function() {
  # n_genes is read off genes
  .assert_error_naming(calc_paralog_norms(matrix(numeric(0), nrow = 2L, ncol = 0L)), "no genes", ERR_EMPTY_INPUT, "genes")
}

# =====================
# filter_paralogs_by_pattern_dosage_effect / _subfunctionalization
# =====================

FILTERS <- list(dosage_effect = filter_paralogs_by_pattern_dosage_effect,
                subfunctionalization = filter_paralogs_by_pattern_subfunctionalization)

test_filter_paralogs_by_pattern <- function() {
  for (name in names(FILTERS)) {
    result <- FILTERS[[name]](.angles(), 0.5, 1L, .one_family(), 1L)
    assert_true(is.integer(result), paste0(name, ": expected integer words"))
    assert_true(identical(dim(result), c(1L, 1L)),
                paste0(name, ": expected (n_mask_chunks, n_families) = (1, 1), got ", toString(dim(result))))
  }
}

test_filter_paralogs_by_pattern_numbers_genes_from_bit_0 <- function() {
  # the one value: angles [0, 1, 0] at 0.5 keep genes 1 and 3 for dosage effect (the word 5),
  # gene 2 for subfunctionalization (the word 2)
  dosage <- filter_paralogs_by_pattern_dosage_effect(.angles(), 0.5, 1L, .one_family(), 1L)
  subfunc <- filter_paralogs_by_pattern_subfunctionalization(.angles(), 0.5, 1L, .one_family(), 1L)
  assert_true(identical(dosage[1L, 1L], 5L) && identical(subfunc[1L, 1L], 2L),
              paste("got", dosage[1L, 1L], "and", subfunc[1L, 1L]))
}

test_filter_paralogs_by_pattern_rejects_more_families_than_genes <- function() {
  for (name in names(FILTERS)) {
    .assert_error_naming(FILTERS[[name]](.angles(), 0.5, 4L, .one_family(), 1L),
                         paste0(name, ": 4 families of 3 genes"), ERR_INVALID_INPUT, "n_families")
  }
}

test_filter_paralogs_by_pattern_rejects_negative_extents <- function() {
  # the binding's own refusal, before it allocates the output by them
  for (name in names(FILTERS)) {
    .assert_error_naming(FILTERS[[name]](.angles(), 0.5, -1L, .one_family(), 1L),
                         paste0(name, ": n_families -1"), NULL, "n_families")
    .assert_error_naming(FILTERS[[name]](.angles(), 0.5, 1L, .one_family(), -1L),
                         paste0(name, ": n_mask_chunks -1"), NULL, "n_mask_chunks")
  }
}

test_filter_paralogs_by_pattern_rejects_too_few_words <- function() {
  .assert_error_naming(filter_paralogs_by_pattern_dosage_effect(.angles(), 0.5, 1L, .one_family(), 0L),
                       "no word for 3 genes", ERR_INVALID_INPUT, "n_mask_chunks")
}

test_filter_paralogs_by_pattern_rejects_a_family_id_of_0 <- function() {
  # these filters take no unassigned genes
  .assert_error_naming(filter_paralogs_by_pattern_subfunctionalization(.angles(), 0.5, 1L, c(1L, 0L, 1L), 1L),
                       "a family id 0", ERR_INVALID_INPUT, "gene_to_fam")
}

test_filter_paralogs_by_pattern_rejects_an_angle_above_pi <- function() {
  .assert_error_naming(filter_paralogs_by_pattern_dosage_effect(c(0, 4, 0), 0.5, 1L, .one_family(), 1L),
                       "an angle of 4", ERR_INVALID_INPUT, "gene_angles")
}

test_filter_paralogs_by_pattern_rejects_a_nan_threshold <- function() {
  .assert_error_naming(filter_paralogs_by_pattern_subfunctionalization(.angles(), NaN, 1L, .one_family(), 1L),
                       "a NaN threshold", ERR_NAN_INF, "threshold")
}

test_filter_paralogs_by_pattern_rejects_no_genes <- function() {
  # n_genes is read off gene_angles
  .assert_error_naming(filter_paralogs_by_pattern_dosage_effect(numeric(0), 0.5, 1L, integer(0), 1L),
                       "no genes", ERR_EMPTY_INPUT, "gene_angles")
}

test_filter_paralogs_by_pattern_rejects_a_mapping_of_other_genes <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(filter_paralogs_by_pattern_dosage_effect(.angles(), 0.5, 1L, c(1L, 1L), 1L),
                       "a mapping of another number of genes", NULL, "gene_to_fam")
}

# =====================
# calc_paralog_subsets_bitmasks_size
# =====================

test_calc_paralog_subsets_bitmasks_size <- function() {
  result <- calc_paralog_subsets_bitmasks_size(3L, 3L, ALL_THREE)
  assert_true(is.list(result) && identical(sort(names(result)), c("max_subset_size", "work_array_size")),
              paste("expected a list of max_subset_size and work_array_size, got", toString(names(result))))
  assert_true(is.integer(result$max_subset_size) && is.integer(result$work_array_size), "expected integers")
}

test_calc_paralog_subsets_bitmasks_size_hands_back_the_cut_size <- function() {
  # the one value: max_subset_size is in and out; 100 comes back as the three kept genes
  result <- calc_paralog_subsets_bitmasks_size(100L, 3L, ALL_THREE)
  assert_true(identical(result$max_subset_size, 3L), paste("got", result$max_subset_size))
}

test_calc_paralog_subsets_bitmasks_size_rejects_a_negative_max_subset_size <- function() {
  .assert_error_naming(calc_paralog_subsets_bitmasks_size(-1L, 3L, ALL_THREE),
                       "max_subset_size -1", ERR_INVALID_INPUT, "max_subset_size")
}

test_calc_paralog_subsets_bitmasks_size_rejects_more_genes_than_bits <- function() {
  .assert_error_naming(calc_paralog_subsets_bitmasks_size(3L, 33L, ALL_THREE),
                       "33 genes in one word", ERR_INVALID_INPUT, "n_genes")
}

test_calc_paralog_subsets_bitmasks_size_rejects_no_genes <- function() {
  .assert_error_naming(calc_paralog_subsets_bitmasks_size(3L, 0L, ALL_THREE),
                       "n_genes 0", ERR_INVALID_INPUT, "n_genes")
}

run_all_tests()
