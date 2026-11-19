# Comprehensive R test suite for tox_paralog_analysis (mirrors python/test/mod_test_tox_paralog_analysis.py)

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

normalize_unit <- function(v) {
  n <- sqrt(sum(v * v))
  if (n == 0) return(v)
  v / n
}

test_paralog_functions <- function() {

  # Testing mask logic
  n_paralogs <- 5L
  i_paralog <- 2L

  chunk_count <- bitmask_n_words(n_paralogs)
  assert_true(chunk_count == 1L)

  bitmask <- integer(chunk_count)
  bitmask[(i_paralog %/% 32L) + 1L] <- bitwShiftL(1L, i_paralog %% 32L)
  state <- bitmask_test(n_paralogs, bitmask, i_paralog + 1L)
  assert_true(isTRUE(state))

  # Testing pattern filtering
  n_families <- 1L
  gene_to_fam <- rep(1L, n_paralogs)
  angles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  threshold <- 0.6

  dosage_mask <- filter_paralogs_by_pattern_dosage_effect(angles, threshold, n_families, gene_to_fam, chunk_count)
  assert_true(all(dim(dosage_mask) == c(1L, 1L)))
  assert_true(dosage_mask[1, 1] == 7L)

  subfunc_mask <- filter_paralogs_by_pattern_subfunctionalization(angles, threshold, n_families, gene_to_fam, chunk_count)
  assert_true(all(dim(subfunc_mask) == c(1L, 1L)))
  assert_true(subfunc_mask[1, 1] == 24L)

  # Testing work-array size calculation
  max_subset_size <- 3L
  work_info <- calc_paralog_subsets_bitmasks_size(max_subset_size, n_paralogs, as.integer(dosage_mask[, 1]))
  assert_true(work_info$work_array_size == 3L)
  assert_true(work_info$max_subset_size == 3L)

  # Testing dosage-effect detection
  ancestor <- c(1.0, 1.0)
  paralogs <- matrix(c(
    1.1, 1.2, 1.3, 1.4, 1.5,
    0.9, 0.8, 0.7, 0.6, 0.5
  ), nrow = 2, byrow = TRUE)

  dosage_result <- detect_dosage_effect(
    ancestor = ancestor,
    genes = paralogs,
    filtered_paralogs_bitmask = as.integer(dosage_mask[, 1]),
    max_subset_size = work_info$max_subset_size,
    gain_gamma = 0.1,
    max_angle = pi
  )

  # the work array is returned whole; n_results says how much of it is filled
  assert_true(dosage_result$n_results == 3L)
  dosage_found <- dosage_result$paralog_subsets_bitmasks[1, seq_len(dosage_result$n_results)]
  assert_true(all(dosage_found == c(6L, 3L, 5L)))

  # an over-large max_subset_size is capped by calc_paralog_subsets_bitmasks_size (inout
  # feedback): the same result comes back at once, not after a runaway subset-size loop
  capped <- detect_dosage_effect(
    ancestor = ancestor,
    genes = paralogs,
    filtered_paralogs_bitmask = as.integer(dosage_mask[, 1]),
    max_subset_size = 2000000000L,
    gain_gamma = 0.1,
    max_angle = pi
  )
  assert_true(capped$n_results == 3L)
  assert_true(all(capped$paralog_subsets_bitmasks[1, seq_len(capped$n_results)] == c(6L, 3L, 5L)))

  # Testing subfunctionalization detection
  norms <- sqrt(colSums(paralogs ^ 2))
  sorted_perm <- as.integer(order(norms))

  subfunc_info <- calc_paralog_subsets_bitmasks_size(max_subset_size, n_paralogs,
                                                     as.integer(subfunc_mask[, 1]))

  subfunc_result <- detect_subfunctionalization(
    ancestor = ancestor,
    genes = paralogs,
    rdi_threshold = 0.5,
    filtered_paralogs_bitmask = as.integer(subfunc_mask[, 1]),
    max_subset_size = subfunc_info$max_subset_size,
    paralog_norms = norms,
    sorted_paralog_norms_perm = sorted_perm
  )

  assert_true(subfunc_result$n_results == 0L)

  # Edge cases. Zero genes need zero words, so bitmask_n_words returns rather than
  # raising; the Python suite asserts no error too.
  assert_true(bitmask_n_words(0L) == 0L)

  single_mask <- bitmask_n_words(1L)
  assert_true(single_mask == 1L)

}


test_a_gene_alone_in_bit_31_survives_the_chain <- function() {
  # Gene 32 is bit 31 of word 1. Kept alone in that word, the word is 0x80000000, which R
  # prints as NA_integer_ (issue #206): every step of filter -> calc -> detect takes it back.
  n_genes <- 33L
  angles <- rep(1.0, n_genes)
  angles[c(32L, 33L)] <- 0.1
  gene_to_fam <- rep(1L, n_genes)

  family_bitmasks <- filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1L, gene_to_fam,
                                                              bitmask_n_words(n_genes))
  kept <- family_bitmasks[, 1]
  assert_true(identical(kept, c(NA_integer_, 1L)), "genes 32 and 33: the NA word, then bit 0 of word 2")

  sizing <- calc_paralog_subsets_bitmasks_size(2L, n_genes, kept)
  assert_true(sizing$max_subset_size == 2L, "two genes are kept, so a subset holds two at most")

  # only the two kept genes together exceed the ancestor's magnitude by 10%
  genes <- matrix(rep(c(0.6, 0.6), n_genes), nrow = 2L)
  found <- detect_dosage_effect(c(1.0, 1.0), genes, kept, sizing$max_subset_size, gain_gamma = 0.1)
  assert_true(found$n_results == 1L, "the pair {32, 33} is the one dosage subset")
  assert_true(identical(found$paralog_subsets_bitmasks[, 1], c(NA_integer_, 1L)),
              "the subset comes back as the same words, the NA word included")

  # an NA in an ordinary integer argument is still refused: it is no family index
  assert_true(inherits(tryCatch(filter_paralogs_by_pattern_dosage_effect(
                angles, 0.5, 1L, c(NA_integer_, gene_to_fam[-1]), bitmask_n_words(n_genes)),
                error = function(e) e), "tox_na_error"),
              "gene_to_fam = NA is refused as NA")
}

test_detect_neofunctionalization <- function() {

  # Case 1: differences below threshold -> all FALSE
  ancestors <- matrix(c(5, 3,
                        2, 1), nrow = 2, ncol = 2)
  ancestors <- apply(ancestors, 2, normalize_unit)

  gene_to_fam <- as.integer(c(1, 2, 1))
  thresholds <- c(0.05, 0.05)

  genes <- matrix(0.0, nrow = 2, ncol = 3)
  for (i in seq_len(3)) {
    genes[, i] <- ancestors[, gene_to_fam[i]]
  }

  neofunc <- detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
  expected <- matrix(FALSE, nrow = 3, ncol = 2)
  assert_true(identical(neofunc, expected))

  # Case 2: differences above threshold -> some TRUE
  thresholds <- c(0.2, 0.2)
  genes2 <- matrix(0.0, nrow = 2, ncol = 3)
  for (i in seq_len(3)) {
    genes2[, i] <- ancestors[, gene_to_fam[i]] - thresholds * gene_to_fam[i]
  }

  neofunc2 <- detect_neofunctionalization(ancestors, genes2, gene_to_fam, thresholds)
  expected2 <- matrix(c(FALSE, TRUE, FALSE,
                        FALSE, TRUE, FALSE), nrow = 3, ncol = 2)
  assert_true(identical(neofunc2, expected2))

}

run_all_tests()
