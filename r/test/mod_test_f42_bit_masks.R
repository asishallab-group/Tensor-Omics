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

test_bitmask_n_words <- function() {
  n_words <- bitmask_n_words(4L)
  assert_true(is.integer(n_words) && length(n_words) == 1L, "expected an integer scalar")
  assert_equal_int(n_words, 1L, "4 flags take 1 word")
  assert_equal_int(bitmask_n_words(32L), 1L, "32 flags fill exactly 1 word")
  assert_equal_int(bitmask_n_words(33L), 2L, "33 flags take 2 words")

  # an empty mask takes no words, as the paralog code's mask_chunk_count did
  assert_equal_int(bitmask_n_words(0L), 0L, "0 flags take 0 words")

  assert_error(bitmask_n_words(-1L), "n_bits = -1 must be rejected", ERR_INVALID_INPUT)
}

test_bitmask_from_logical <- function() {
  bitmask <- bitmask_from_logical(1L, FLAGS_1011)
  assert_true(is.integer(bitmask), "expected an integer vector")
  assert_equal_int(bitmask, WORD_1011, "[T, T, F, T] must pack to 11")

  # flag 33 is bit 0 of the second word
  flags <- logical(33L)
  flags[33L] <- TRUE
  assert_equal_int(bitmask_from_logical(2L, flags), c(0L, 1L), "flag 33 must be bit 0 of word 2")

  assert_error(bitmask_from_logical(2L, FLAGS_1011), "n_words = 2 for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bitmask_from_logical(0L, FLAGS_1011), "n_words = 0 for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bitmask_from_logical(0L, logical(0L)), "no flags (n_bits = 0) must be rejected",
               ERR_INVALID_INPUT)
}

test_bitmask_to_logical <- function() {
  flags <- bitmask_to_logical(4L, WORD_1011)
  assert_true(is.logical(flags), "expected a logical vector")
  assert_true(identical(flags, FLAGS_1011), "11 must unpack to [T, T, F, T]")

  assert_error(bitmask_to_logical(0L, 0L), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bitmask_to_logical(4L, c(WORD_1011, 0L)), "2 words for 4 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bitmask_to_logical(33L, WORD_1011), "1 word for 33 flags must be rejected",
               ERR_INVALID_INPUT)
}

test_bitmask_test <- function() {
  is_set <- bitmask_test(4L, WORD_1011, 4L)
  assert_true(is.logical(is_set) && length(is_set) == 1L, "expected a logical scalar")
  assert_true(is_set, "flag 4 of 11 is set")
  assert_false(bitmask_test(4L, WORD_1011, 3L), "flag 3 of 11 is clear")

  assert_error(bitmask_test(4L, WORD_1011, 0L), "i_bit = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bitmask_test(4L, WORD_1011, 5L), "i_bit = n_bits + 1 must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bitmask_test(0L, WORD_1011, 1L), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bitmask_test(4L, c(WORD_1011, 0L), 1L), "2 words for 4 flags must be rejected",
               ERR_INVALID_INPUT)
}

test_bitmasks_from_logical_2D <- function() {
  # 33 flags per mask take 2 words, and every column starts on a word of its own: column 2's
  # flag 1 is bit 0 of its own first word, not bit 1 of column 1's second word
  flags <- matrix(FALSE, nrow = 33L, ncol = 2L)
  flags[, 1L] <- TRUE
  flags[c(1L, 33L), 2L] <- TRUE
  bitmasks <- bitmasks_from_logical_2D(2L, flags)
  assert_true(is.matrix(bitmasks) && is.integer(bitmasks), "expected an integer matrix")
  assert_equal_int(dim(bitmasks), c(2L, 2L), "expected dim (2, 2)")
  # column 1: 32 set flags are all bits of word 1 (-1 as int32), then flag 33 alone
  assert_equal_int(bitmasks[, 1L], c(-1L, 1L), "column 1 must be [-1, 1]")
  assert_equal_int(bitmasks[, 2L], c(1L, 1L), "column 2 must be [1, 1]")

  assert_error(bitmasks_from_logical_2D(1L, flags), "n_words = 1 for 33 flags must be rejected",
               ERR_INVALID_INPUT)
  assert_error(bitmasks_from_logical_2D(1L, matrix(logical(0L), nrow = 0L, ncol = 2L)),
               "no flags (n_bits = 0) must be rejected", ERR_INVALID_INPUT)
}

test_bitmasks_to_logical_2D <- function() {
  # two masks of 4 flags: one word each
  bitmasks <- matrix(c(WORD_1011, 0L), nrow = 1L)
  flags <- bitmasks_to_logical_2D(4L, bitmasks)
  assert_true(is.matrix(flags) && is.logical(flags), "expected a logical matrix")
  assert_equal_int(dim(flags), c(4L, 2L), "expected dim (4, 2)")
  assert_true(identical(flags[, 1L], FLAGS_1011), "column 1 must unpack to [T, T, F, T]")
  assert_false(any(flags[, 2L]), "column 2 is the empty mask")

  # the round trip through the padded layout
  padded_flags <- matrix(FALSE, nrow = 33L, ncol = 2L)
  padded_flags[seq(1L, 31L, by = 3L), 1L] <- TRUE
  padded_flags[33L, 2L] <- TRUE
  round_trip <- bitmasks_to_logical_2D(33L, bitmasks_from_logical_2D(2L, padded_flags))
  assert_true(identical(round_trip, padded_flags), "unpacking must undo packing")

  assert_error(bitmasks_to_logical_2D(0L, bitmasks), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
  assert_error(bitmasks_to_logical_2D(33L, bitmasks),
               "1 word per mask for 33 flags must be rejected", ERR_INVALID_INPUT)
}

test_layout_matches_int_to_bits <- function() {
  # the module documentation promises intToBits decodes the words directly. Flags 32 and 64 make
  # both full words negative, which R handles; only a word with bit 31 alone is NA (issue #206)
  n_bits <- 70L
  flags <- logical(n_bits)
  flags[c(1L, 6L, 32L, 33L, 41L, 64L, 65L, 70L)] <- TRUE
  bitmask <- bitmask_from_logical(bitmask_n_words(n_bits), flags)
  assert_true(identical(decode_with_int_to_bits(bitmask, n_bits), bitmask_to_logical(n_bits, bitmask)),
              "intToBits must decode what bitmask_to_logical decodes")
  assert_true(identical(decode_with_int_to_bits(bitmask, n_bits), flags), "intToBits must recover the flags")
}

test_word_with_only_bit_31_set_round_trips <- function() {
  # Flag 32 alone is the word 0x80000000, the most negative int32, and R uses exactly that bit
  # pattern for NA_integer_ (issue #206). It prints as NA, is still that word, and passes back.
  flags <- logical(32L)
  flags[32L] <- TRUE
  bitmask <- bitmask_from_logical(1L, flags)
  assert_true(identical(bitmask, NA_integer_), "flag 32 alone is the word R prints as NA")
  assert_true(identical(decode_with_int_to_bits(bitmask, 32L), flags),
              "intToBits decodes the NA word as flag 32 alone")

  assert_true(identical(bitmask_test(32L, bitmask, 32L), TRUE), "flag 32 of the NA word is set")
  assert_true(identical(bitmask_test(32L, bitmask, 31L), FALSE), "flag 31 of the NA word is not")
  assert_true(identical(bitmask_to_logical(32L, bitmask), flags), "the NA word unpacks to flag 32 alone")

  # beside other words, and in the second word of a mask
  flags_40 <- logical(40L)
  flags_40[c(32L, 33L)] <- TRUE
  words <- bitmask_from_logical(2L, flags_40)
  assert_true(identical(words, c(NA_integer_, 1L)), "flag 32 alone in word 1, flag 33 in word 2")
  assert_true(identical(bitmask_to_logical(40L, words), flags_40), "a mask holding the NA word round-trips")
}

test_word_with_only_bit_31_set_round_trips_in_2D <- function() {
  # the 2-D form takes the same word in any column
  flags <- matrix(FALSE, nrow = 32L, ncol = 3L)
  flags[32L, 1L] <- TRUE
  flags[c(1L, 2L, 4L), 2L] <- TRUE
  flags[32L, 3L] <- TRUE
  bitmasks <- bitmasks_from_logical_2D(1L, flags)
  assert_true(identical(bitmasks, matrix(c(NA_integer_, WORD_1011, NA_integer_), nrow = 1L)),
              "columns 1 and 3 are the NA word, column 2 is 11")
  assert_true(identical(bitmasks_to_logical_2D(32L, bitmasks), flags), "every column round-trips")
}

test_a_double_NA_is_no_word <- function() {
  # a double has a missing value of its own, which is no int32 bit pattern, so it is refused,
  # with a hint how to pass the word with only bit 31 set
  err <- tryCatch(bitmask_to_logical(32L, NA_real_), error = function(e) e)
  assert_true(inherits(err, "tox_na_error"), "a double NA is refused as NA")
  assert_true(grepl("integer vector", conditionMessage(err), fixed = TRUE),
              "the refusal says to pass the words as an integer vector")
  assert_error(bitmasks_to_logical_2D(32L, matrix(c(1, NA), nrow = 1L)), "a double NA in a column is refused")

  # whole doubles in the integer range are words, as for any integer array
  assert_true(identical(bitmask_to_logical(4L, 11), FLAGS_1011), "the double 11 is the word 11")
}

run_all_tests()
