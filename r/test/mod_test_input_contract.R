# The input contract of the R binding: every case of the table, accepted and refused.
#
# The table is the one both bindings follow (helper/codegen/design/language-layers.md, "what an
# argument accepts"); python/test/mod_test_input_contract.py is the same list in Python. A
# refusal is the binding's own: a tox_input_error condition naming the argument, code NA (#221).
# Its elements are checked in the C shim as it converts them, its kind and rank in R.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

EXPR <- matrix(c(1, 2, 3, 4), nrow = 1)          # one axis, four genes
FAMILY <- c(1L, 1L, 1L, 1L)
MASK <- c(TRUE, FALSE, TRUE, FALSE)
SCRATCH <- tempfile("tox_contract_")
dir.create(SCRATCH)
scratch <- function(name) file.path(SCRATCH, name)

# `expr` raises the binding's own refusal of `argument`, as `class`, its message matching `pattern`
refused <- function(expr, argument, pattern, class = "tox_type_error") {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(sprintf("'%s' was not refused", argument), call. = FALSE)
  message <- conditionMessage(err)
  if (!inherits(err, class) || !inherits(err, "tox_input_error"))
    stop(sprintf("expected a %s, got <%s>: %s", class, paste(class(err), collapse = ", "), message),
         call. = FALSE)
  if (!is.na(err$code)) stop(sprintf("a refusal carries no code, got %d", err$code), call. = FALSE)
  if (!identical(err$argument, argument))
    stop(sprintf("expected '%s' named, got '%s': %s", argument, err$argument, message), call. = FALSE)
  if (!grepl(pattern, message)) stop(sprintf("'%s' does not match '%s'", message, pattern), call. = FALSE)
  invisible(TRUE)
}

centroid <- function(gene_to_family = FAMILY, n_families = 1L, mask = MASK, expr = EXPR)
  group_centroid_orthologs(expr, gene_to_family, n_families, mask)

# =====================================================================================
# integer(int32)
# =====================================================================================

test_integer_accepts_integers_and_whole_doubles <- function() {
  # R's own 2 is a double: 96 calls in the suites pass one (#221 A6.4)
  assert_equal_numeric(centroid(c(1, 1, 1, 1), n_families = 1), 2)
  assert_equal_numeric(centroid(FAMILY, n_families = 1L), 2)
  # -0 is 0: a gene in no family, which the mask leaves out anyway
  assert_equal_numeric(centroid(c(1, -0, 1, 1)), 2)
  # a compact sequence, which R knows holds no NA
  assert_true(is.null(validate_gene_to_family_mapping(1:4, 4L)))
}

test_integer_refuses_logicals <- function() {
  refused(centroid(n_families = TRUE), "n_families", "not logical")
  refused(centroid(c(TRUE, TRUE, TRUE, TRUE)), "gene_to_family", "not logical")
}

test_integer_refuses_fractions_nan_na_and_infinities <- function() {
  # 1.7 was truncated to 1 before (#221 item 3)
  refused(centroid(c(1.7, 1, 1, 1)), "gene_to_family", "element 1 is 1.7")
  refused(centroid(c(1, NaN, 1, 1)), "gene_to_family", "element 2 is NaN")
  refused(centroid(c(1, 1, Inf, 1)), "gene_to_family", "element 3 is Inf")
  refused(centroid(c(1, NA, 1, 1)), "gene_to_family", "contains NA", "tox_na_error")
  refused(centroid(c(1L, NA, 1L, 1L)), "gene_to_family", "contains NA", "tox_na_error")
  refused(centroid(n_families = 1.5), "n_families", "not 1.5")
  refused(centroid(n_families = NaN), "n_families", "not NaN")
  refused(centroid(n_families = NA_integer_), "n_families", "contains NA", "tox_na_error")
}

test_integer_refuses_out_of_range <- function() {
  refused(centroid(c(2^31, 1, 1, 1)), "gene_to_family", "element 1 is 2147483648")
  refused(centroid(c(1, -2^31, 1, 1)), "gene_to_family", "element 2 is -2147483648")
  refused(centroid(n_families = 3e9), "n_families", "not 3e\\+09")
}

test_integer_refuses_strings_complex_factor_date_and_integer64 <- function() {
  refused(centroid(c("1", "1", "1", "1")), "gene_to_family", "not character")
  refused(centroid(complex(real = FAMILY)), "gene_to_family", "not complex")
  # as.integer() would give a factor's codes, not its labels
  refused(centroid(factor(c(1, 1, 1, 1))), "gene_to_family", "not factor")
  refused(centroid(Sys.Date() + 0:3), "gene_to_family", "not Date")
  # bit64's storage, without needing bit64: doubles holding int64 bits
  int64 <- structure(c(4.9e-324, 4.9e-324, 4.9e-324, 4.9e-324), class = "integer64")
  refused(centroid(int64), "gene_to_family", "integer64.*as.integer\\(\\)")
}

# =====================================================================================
# real(real64)
# =====================================================================================

test_real_accepts_integers_and_doubles <- function() {
  assert_equal_numeric(euclidean_distance(c(3L, 4L), c(0, 0)), 5)
  assert_equal_numeric(euclidean_distance(c(3, 4), c(0L, 0L)), 5)
  assert_equal_numeric(calc_quantile(c(1, 2, 3), 1L), 3)
}

test_real_leaves_nan_to_the_procedure <- function() {
  # DM_ALLOW_NAN decides, so the binding passes it on and the Fortran refuses it with a code
  err <- tryCatch({ euclidean_distance(c(NaN, 4), c(0, 0)); NULL }, error = function(e) e)
  assert_true(!is.null(err) && !is.na(err$code), "NaN reaches the procedure")
}

test_real_refuses_logicals_strings_complex_and_dates <- function() {
  refused(euclidean_distance(c(TRUE, FALSE), c(0, 0)), "vec1", "not logical")
  refused(euclidean_distance(c("1.5", "2"), c(0, 0)), "vec1", "not character")
  refused(euclidean_distance(complex(real = c(1, 2)), c(0, 0)), "vec1", "not complex")
  refused(euclidean_distance(Sys.Date() + 0:1, c(0, 0)), "vec1", "not Date")
  refused(euclidean_distance(Sys.time() + 0:1, c(0, 0)), "vec1", "not POSIXct")
  refused(euclidean_distance(list(3, 4), c(0, 0)), "vec1", "not list")
  refused(calc_quantile(c(1, 2), TRUE), "level", "not logical")
  refused(calc_quantile(c(1, 2), "0.5"), "level", "not character")
}

# =====================================================================================
# logical
# =====================================================================================

test_logical_accepts_logicals_and_numbers_that_are_0_or_1 <- function() {
  # c(1, 0, 1, 0) is R's own literal for a 0/1 mask, and Python accepts its [1, 0, 1, 0]
  for (mask in list(MASK, c(1L, 0L, 1L, 0L), c(1, 0, 1, 0), c(1, -0, 1, 0)))
    assert_equal_numeric(centroid(mask = mask), 2)
  negative <- matrix(c(-1, 1), nrow = 1)
  for (flag in list(TRUE, 1L, 1))
    assert_error(validate_expression_data(negative, flag), "checked", ERR_INVALID_INPUT)
  for (flag in list(FALSE, 0L, 0))
    assert_true(is.null(validate_expression_data(negative, flag)))
}

test_logical_refuses_other_numbers_na_and_strings <- function() {
  refused(centroid(mask = c(1L, 2L, 3L, 4L)), "ortholog_set", "element 2 is 2")
  refused(centroid(mask = c(0.5, 0, 0, 0)), "ortholog_set", "element 1 is 0.5")
  refused(centroid(mask = c(1, NaN, 0, 0)), "ortholog_set", "element 2 is NaN")
  refused(centroid(mask = c(TRUE, NA, TRUE, FALSE)), "ortholog_set", "contains NA", "tox_na_error")
  refused(centroid(mask = c(1L, NA, 1L, 0L)), "ortholog_set", "contains NA", "tox_na_error")
  refused(centroid(mask = c("no", "", "", "")), "ortholog_set", "not character")
  refused(centroid(mask = list(TRUE, FALSE, TRUE, FALSE)), "ortholog_set", "not list")
  refused(validate_expression_data(EXPR, 2), "check_non_negative", "not 2")
  refused(validate_expression_data(EXPR, NA), "check_non_negative", "contains NA", "tox_na_error")
  refused(validate_expression_data(EXPR, "FALSE"), "check_non_negative", "not character")
}

test_a_serialized_logical_array_refuses_other_numbers <- function() {
  # c(2L, 3L) and c(0.5, -1) came back TRUE TRUE before (#221 item 5)
  refused(serialize_logical_helper(c(2L, 3L), scratch("l.bin")), "arr", "element 1 is 2")
  refused(serialize_logical_helper(c(0.5, -1), scratch("l.bin")), "arr", "element 1 is 0.5")
  serialize_logical_helper(c(1, 0, 1), scratch("l.bin"))
  assert_true(identical(as.vector(deserialize_logical_helper(scratch("l.bin"))), c(TRUE, FALSE, TRUE)))
}

# =====================================================================================
# a string
# =====================================================================================

test_a_string_accepts_one_string <- function() {
  serialize_char_helper("x", scratch("s.bin"))
  assert_true(identical(as.vector(deserialize_char_helper(scratch("s.bin"))), "x"))
}

test_a_string_refuses_a_vector_na_numbers_and_a_factor <- function() {
  # only the first file was written before (#221 N3, G-S5)
  refused(serialize_char_helper("x", c(scratch("f1.bin"), scratch("f2.bin"))), "filename", "single string")
  assert_false(file.exists(scratch("f1.bin")))
  refused(serialize_char_helper("x", character(0)), "filename", "single string")
  refused(serialize_char_helper("x", NA_character_), "filename", "contains NA", "tox_na_error")
  refused(serialize_char_helper("x", 42), "filename", "not numeric")
  refused(serialize_char_helper("x", factor(scratch("fac.bin"))), "filename", "not factor")
}

test_a_string_that_is_not_utf8_is_refused <- function() {
  refused(serialize_char_helper("x", "f\xffg.bin"), "filename", "element 1 is not valid UTF-8")
  # an encoded surrogate is not well-formed UTF-8 either
  refused(validate_string_array_uniqueness(c("a", "\xed\xa0\x80")), "str_arr", "element 2 is not valid UTF-8")
  # an overlong form, a code point past U+10FFFF, a truncated sequence
  for (bad in c("\xc0\xaf", "\xf4\x90\x80\x80", "\xe2\x82"))
    refused(validate_string_array_uniqueness(c(bad, "b")), "str_arr", "element 1 is not valid UTF-8")
  # and what is well-formed passes, up to four bytes
  assert_true(is.null(validate_string_array_uniqueness(c("\u00e9", "\u20ac", "\U0001F600"))))
}

test_a_string_longer_than_its_declared_length_is_refused_not_cut <- function() {
  # `delimiter = "::"` reached Fortran as ":" before
  refused(read_expression_vectors_tsv(scratch("none.tsv"), "g", matrix(0, 1, 1), 0L, 1L, 2L, 1L,
                                      delimiter = "::"),
          "delimiter", "element 1 is 2 bytes long in UTF-8, but it is declared 1 long")
}

# =====================================================================================
# an array of strings
# =====================================================================================

test_strings_accept_a_character_vector_and_a_factor_by_its_labels <- function() {
  assert_true(is.null(validate_string_array_uniqueness(c("a", "b"))))
  assert_true(is.null(validate_string_array_uniqueness(factor(c("a", "b")))))
  # its labels, not its codes: two equal labels are duplicates
  assert_error(validate_string_array_uniqueness(factor(c("a", "a"))), "labels", ERR_INVALID_INPUT)
}

test_strings_refuse_na_numbers_and_a_matrix <- function() {
  refused(validate_string_array_uniqueness(c("a", NA)), "str_arr", "contains NA", "tox_na_error")
  refused(validate_string_array_uniqueness(c(1, 2)), "str_arr", "not numeric")
  refused(validate_string_array_uniqueness(matrix(c("a", "b", "c", "d"), 2)), "str_arr", "not matrix")
}

# =====================================================================================
# encodings
# =====================================================================================

test_a_latin1_string_reaches_fortran_as_utf8 <- function() {
  utf8 <- "\u00e9"
  latin1 <- iconv(utf8, "UTF-8", "latin1")
  assert_true(Encoding(latin1) == "latin1", "the probe's premise")
  # the same string to Fortran: a duplicate (#221 item 2)
  assert_error(validate_string_array_uniqueness(c(utf8, latin1)), "one string", ERR_INVALID_INPUT)
  serialize_char_helper(latin1, scratch("latin1.bin"))
  back <- deserialize_char_helper(scratch("latin1.bin"))
  assert_true(identical(charToRaw(back), as.raw(c(0xc3, 0xa9))), "UTF-8 bytes")
}

test_a_returned_string_is_marked_utf8 <- function() {
  # it read as unknown before, which a non-UTF-8 locale takes as native (#221 N6)
  serialize_char_helper(c("\u00e9t\u00e9", "plain"), scratch("utf8.bin"))
  back <- as.vector(deserialize_char_helper(scratch("utf8.bin")))
  assert_true(identical(back, c("\u00e9t\u00e9", "plain")))
  assert_true(Encoding(back)[1] == "UTF-8", "marked UTF-8")
}

# What R marks a string as decides how it crosses, never the session's locale: in the C locale an
# unmarked "été" was translated into the escape text "<c3><a9>t<c3><a9>", and that passed
# (#221 WP8). Run in the C locale and in the session's own, a UTF-8 one wherever the suites run.
in_each_locale <- function(check) {
  original <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", original))
  for (locale in unique(c("C", original))) {
    if (Sys.setlocale("LC_CTYPE", locale) == "") stop(sprintf("no locale %s", locale), call. = FALSE)
    check(locale)
  }
}

test_a_string_reaches_fortran_by_its_mark_in_any_locale <- function() {
  unmarked <- rawToChar(as.raw(c(0xc3, 0xa9, 0x74, 0xc3, 0xa9)))  # UTF-8 bytes, as readLines gives them
  marked <- unmarked
  Encoding(marked) <- "UTF-8"
  latin1 <- iconv(marked, "UTF-8", "latin1")
  in_each_locale(function(locale) {
    for (spelling in list(unmarked, marked, latin1)) {
      serialize_char_helper(spelling, scratch("mark.bin"))
      back <- deserialize_char_helper(scratch("mark.bin"))
      assert_true(identical(charToRaw(back), charToRaw(marked)),
                  sprintf("%s in %s: %s", Encoding(spelling), locale, paste(charToRaw(back), collapse = " ")))
    }
    # one string to Fortran, however it is spelled: duplicates
    assert_error(validate_string_array_uniqueness(c(unmarked, marked, latin1)), locale, ERR_INVALID_INPUT)
  })
}

test_a_string_that_cannot_cross_as_utf8_is_refused_in_any_locale <- function() {
  native <- rawToChar(as.raw(c(0x66, 0xe9, 0x67)))                 # latin1 bytes, unmarked
  bytes <- "é"
  Encoding(bytes) <- "bytes"
  undefined <- rawToChar(as.raw(c(0x61, 0x81)))                    # no character in Windows-1252
  Encoding(undefined) <- "latin1"
  in_each_locale(function(locale) {
    refused(serialize_char_helper(native, scratch("native.bin")), "arr", "element 1 is not valid UTF-8")
    refused(serialize_char_helper(bytes, scratch("bytes.bin")), "arr", "element 1 is marked as bytes")
    refused(serialize_char_helper(undefined, scratch("undefined.bin")), "arr",
            "element 1 is marked latin1, but holds a byte that is no character in it")
  })
}

# =====================================================================================
# trailing blanks: Fortran semantics, accepted and dropped
# =====================================================================================

test_trailing_blanks_are_accepted_and_not_significant <- function() {
  assert_error(validate_string_array_uniqueness(c("a", "a ")), "duplicates", ERR_INVALID_INPUT)
  serialize_char_helper(c("ab  ", " c"), scratch("blanks.bin"))
  assert_true(identical(as.vector(deserialize_char_helper(scratch("blanks.bin"))), c("ab", " c")))
}

# =====================================================================================
# modes
# =====================================================================================

DISTANCES <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3)

test_a_mode_matches_whole_with_the_case_rule_kept <- function() {
  assert_true(length(linkage_clustering(DISTANCES, "ward")$heights) == 2L)
  assert_true(length(linkage_clustering(DISTANCES, "WARD")$heights) == 2L)
  # trailing blanks are Fortran semantics, as for every string; a leading one is part of it
  assert_true(length(linkage_clustering(DISTANCES, "ward  ")$heights) == 2L)
  refused(linkage_clustering(DISTANCES, " ward"), "method", "not \" ward\"")
}

test_a_mode_is_never_truncated_and_the_message_lists_the_modes <- function() {
  refused(linkage_clustering(DISTANCES, "weightedXYZ_not_a_method"), "method",
          "must be one of \"average\", .*\"weighted\".*not \"weightedxyz_not_a_method\"")
  refused(linkage_clustering(DISTANCES, 5), "method", "single string")
}

# =====================================================================================
# rank
# =====================================================================================

test_rank_must_be_the_declared_one_and_nothing_is_flattened <- function() {
  # a matrix for a vector was flattened before (#221 G-I3)
  refused(euclidean_distance(matrix(c(3, 4), 1), c(0, 0)), "vec1", "numeric vector, not matrix")
  refused(centroid(matrix(FAMILY, 2)), "gene_to_family", "integer vector, not matrix")
  refused(centroid(mask = matrix(MASK, 2)), "ortholog_set", "logical vector, not matrix")
  refused(centroid(expr = c(1, 2, 3, 4)), "expression_vectors", "numeric matrix")
  refused(centroid(n_families = c(1L, 1L)), "n_families", "integer scalar")
}

test_a_one_dimensional_array_or_table_is_a_vector <- function() {
  assert_equal_numeric(euclidean_distance(array(c(3, 4)), c(0, 0)), 5)
  assert_equal_numeric(centroid(array(FAMILY)), 2)
  counts <- table(c("a", "a", "a", "b", "b", "b", "b"))       # 3 and 4
  assert_equal_numeric(euclidean_distance(counts, c(0, 0)), 5)
}

test_a_numeric_data_frame_is_a_matrix <- function() {
  assert_equal_numeric(centroid(expr = as.data.frame(EXPR)), 2)
  refused(centroid(expr = data.frame(a = 1, b = "x")), "expression_vectors", "numeric matrix, not data.frame")
}

test_a_serialized_integer_array_takes_any_rank_and_checks_its_values <- function() {
  serialize_int_helper(matrix(c(1, 2, 3, 4), 2), scratch("ints2d.bin"))
  back <- deserialize_int_helper(scratch("ints2d.bin"))
  assert_true(identical(back, matrix(1:4, 2)))
  refused(serialize_int_helper(c(1.5, 2), scratch("frac.bin")), "arr", "element 1 is 1.5")
}

run_all_tests()
unlink(SCRATCH, recursive = TRUE)
