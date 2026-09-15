source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# The values are tested in Fortran. These tests check that each procedure can be called, its
# return type and shape, and its documented errors, plus one hand-derived value per procedure:
# it proves the binding passes the layout through unchanged. Flag i is bit mod(i-1, 32) of word
# (i-1) %/% 32 + 1, least significant first, so the flags [T, T, F, T] are the word 0b1011 = 11.
FLAGS_1011 <- c(TRUE, TRUE, FALSE, TRUE)
WORD_1011 <- 11L

# The published layout, decoded the way the module documentation says R does it.
decode_with_int_to_bits <- function(words, n_bits) as.logical(intToBits(words))[seq_len(n_bits)]

test_bit_mask_n_words <- function() {
  n_words <- bit_mask_n_words(4L)
  assert_true(is.integer(n_words) && length(n_words) == 1L, "expected an integer scalar")
  assert_equal_int(n_words, 1L, "4 flags take 1 word")
  assert_equal_int(bit_mask_n_words(32L), 1L, "32 flags fill exactly 1 word")
  assert_equal_int(bit_mask_n_words(33L), 2L, "33 flags take 2 words")

  # an empty mask takes no words, as the paralog code's mask_chunk_count did
  assert_equal_int(bit_mask_n_words(0L), 0L, "0 flags take 0 words")

  assert_error(bit_mask_n_words(-1L), "n_bits = -1 must be rejected", ERR_INVALID_INPUT)
}

test_bit_mask_from_logical <- function() {
  bit_mask <- bit_mask_from_logical(1L, FLAGS_1011)
  assert_true(is.integer(bit_mask), "expected an integer vector")
  assert_equal_int(bit_mask, WORD_1011, "[T, T, F, T] must pack to 11")

  # flag 33 is bit 0 of the second word
  flags <- logical(33L)
  flags[33L] <- TRUE
  assert_equal_int(bit_mask_from_logical(2L, flags), c(0L, 1L), "flag 33 must be bit 0 of word 2")

  assert_error(bit_mask_from_logical(2L, FLAGS_1011), "n_words = 2 for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bit_mask_from_logical(0L, FLAGS_1011), "n_words = 0 for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bit_mask_from_logical(0L, logical(0L)), "no flags (n_bits = 0) must be rejected",
               ERR_INVALID_INPUT)
}

test_bit_mask_to_logical <- function() {
  flags <- bit_mask_to_logical(4L, WORD_1011)
  assert_true(is.logical(flags), "expected a logical vector")
  assert_true(identical(flags, FLAGS_1011), "11 must unpack to [T, T, F, T]")

  assert_error(bit_mask_to_logical(0L, 0L), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bit_mask_to_logical(4L, c(WORD_1011, 0L)), "2 words for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bit_mask_to_logical(33L, WORD_1011), "1 word for 33 flags must be rejected",
               ERR_INVALID_INPUT)
}

test_bit_mask_test <- function() {
  is_set <- bit_mask_test(4L, WORD_1011, 4L)
  assert_true(is.logical(is_set) && length(is_set) == 1L, "expected a logical scalar")
  assert_true(is_set, "flag 4 of 11 is set")
  assert_false(bit_mask_test(4L, WORD_1011, 3L), "flag 3 of 11 is clear")

  assert_error(bit_mask_test(4L, WORD_1011, 0L), "i_bit = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bit_mask_test(4L, WORD_1011, 5L), "i_bit = n_bits + 1 must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bit_mask_test(0L, WORD_1011, 1L), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bit_mask_test(4L, c(WORD_1011, 0L), 1L), "2 words for 4 flags must be rejected",
               ERR_INVALID_INPUT)
}

test_bit_masks_from_logical_2D <- function() {
  # 33 flags per mask take 2 words, and every column starts on a word of its own: column 2's
  # flag 1 is bit 0 of its own first word, not bit 1 of column 1's second word
  flags <- matrix(FALSE, nrow = 33L, ncol = 2L)
  flags[, 1L] <- TRUE
  flags[c(1L, 33L), 2L] <- TRUE
  bit_masks <- bit_masks_from_logical_2D(2L, flags)
  assert_true(is.matrix(bit_masks) && is.integer(bit_masks), "expected an integer matrix")
  assert_equal_int(dim(bit_masks), c(2L, 2L), "expected dim (2, 2)")
  # column 1: 32 set flags are all bits of word 1 (-1 as int32), then flag 33 alone
  assert_equal_int(bit_masks[, 1L], c(-1L, 1L), "column 1 must be [-1, 1]")
  assert_equal_int(bit_masks[, 2L], c(1L, 1L), "column 2 must be [1, 1]")

  assert_error(bit_masks_from_logical_2D(1L, flags), "n_words = 1 for 33 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bit_masks_from_logical_2D(1L, matrix(logical(0L), nrow = 0L, ncol = 2L)),
               "no flags (n_bits = 0) must be rejected", ERR_INVALID_INPUT)
}

test_bit_masks_to_logical_2D <- function() {
  # two masks of 4 flags: one word each
  bit_masks <- matrix(c(WORD_1011, 0L), nrow = 1L)
  flags <- bit_masks_to_logical_2D(4L, bit_masks)
  assert_true(is.matrix(flags) && is.logical(flags), "expected a logical matrix")
  assert_equal_int(dim(flags), c(4L, 2L), "expected dim (4, 2)")
  assert_true(identical(flags[, 1L], FLAGS_1011), "column 1 must unpack to [T, T, F, T]")
  assert_false(any(flags[, 2L]), "column 2 is the empty mask")

  # the round trip through the padded layout
  padded_flags <- matrix(FALSE, nrow = 33L, ncol = 2L)
  padded_flags[seq(1L, 31L, by = 3L), 1L] <- TRUE
  padded_flags[33L, 2L] <- TRUE
  round_trip <- bit_masks_to_logical_2D(33L, bit_masks_from_logical_2D(2L, padded_flags))
  assert_true(identical(round_trip, padded_flags), "unpacking must undo packing")

  assert_error(bit_masks_to_logical_2D(0L, bit_masks), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bit_masks_to_logical_2D(33L, bit_masks),
               "1 word per mask for 33 flags must be rejected", ERR_INVALID_INPUT)
}

test_layout_matches_int_to_bits <- function() {
  # the module documentation promises intToBits decodes the words directly. Flags 32 and 64 make
  # both full words negative, which R handles; only a word with bit 31 alone is NA (issue #206)
  n_bits <- 70L
  flags <- logical(n_bits)
  flags[c(1L, 6L, 32L, 33L, 41L, 64L, 65L, 70L)] <- TRUE
  bit_mask <- bit_mask_from_logical(bit_mask_n_words(n_bits), flags)
  assert_true(identical(decode_with_int_to_bits(bit_mask, n_bits), bit_mask_to_logical(n_bits, bit_mask)),
              "intToBits must decode what bit_mask_to_logical decodes")
  assert_true(identical(decode_with_int_to_bits(bit_mask, n_bits), flags), "intToBits must recover the flags")
}

test_word_with_only_bit_31_set_is_NA <- function() {
  # KNOWN ISSUE #206. Flag 32 alone is the word 0x80000000, the most negative int32, and R uses
  # exactly that bit pattern for NA_integer_. This test pins today's behaviour, so that fixing
  # #206 makes it fail on purpose: rewrite it then to assert the flag round-trips.
  flags <- logical(32L)
  flags[32L] <- TRUE
  bit_mask <- bit_mask_from_logical(1L, flags)
  assert_true(identical(bit_mask, NA_integer_), "#206: flag 32 alone currently packs to NA_integer_")

  # the bits are all there: intToBits reads NA_integer_'s pattern as flag 32 alone
  assert_true(identical(decode_with_int_to_bits(bit_mask, 32L), flags),
              "#206: intToBits still decodes the NA word as flag 32 alone")

  # but the binding refuses the word on the way back in, before the library is called
  assert_true(inherits(tryCatch(bit_mask_test(32L, bit_mask, 32L), error = function(e) e), "tox_na_error"),
              "#206: bit_mask_test currently rejects the NA word with a tox_na_error")
  assert_error(bit_mask_test(32L, bit_mask, 32L), "#206: bit_mask_test rejects the NA word")
  assert_error(bit_mask_to_logical(32L, bit_mask), "#206: bit_mask_to_logical rejects the NA word")

  # the 2-D form hits the same word in any column
  bit_masks <- bit_masks_from_logical_2D(1L, cbind(flags, FLAGS_1011[c(1:4, rep(3L, 28L))]))
  assert_true(is.na(bit_masks[1L, 1L]) && identical(bit_masks[1L, 2L], WORD_1011),
              "#206: only the column holding 0x80000000 is NA")
  assert_error(bit_masks_to_logical_2D(32L, bit_masks), "#206: bit_masks_to_logical_2D rejects the NA word")
}

run_all_tests()
