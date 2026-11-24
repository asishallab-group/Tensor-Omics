# Generated. Do not edit.

#' The number of words a mask of `n_bits` flags takes.
#'
#' Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmask_n_words}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags in the mask; no flags take no words
#'   The minimum valid value is `0`.
#' @return a integer scalar. number of words the mask takes
#' @export
bitmask_n_words <- function(n_bits) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    .result <- .Call("bitmask_n_words_call", n_bits)
    .arguments <- c("n_bits", "n_words", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_words
}

#' Whether flag `i_bit` of a mask is set.
#'
#' Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmask_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags in the mask
#'   The minimum valid value is `1`.
#' @param bitmask a integer vector. the mask
#' @param i_bit a integer scalar. index of the flag, counting from 1
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_bits`.
#' @return a logical scalar. `TRUE` if the flag is set
#' @export
bitmask_test <- function(n_bits, bitmask, i_bit) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bitmask <- .tox_as_integer_vector(bitmask, "bitmask")
    i_bit <- .tox_as_integer_scalar(i_bit, "i_bit")
    .result <- .Call("bitmask_test_call", n_bits, bitmask, i_bit)
    .arguments <- c("n_bits", "n_words", "bitmask", "i_bit", "is_set", "ierr")
    .sources <- c(NA_character_, "bitmask", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_set
}

#' Packs a `logical(c_bool)` array into a bit mask.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmask_from_logical}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_words a integer scalar. number of words of the mask
#'   The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#'   The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#' @param flags a logical vector. the flags to pack
#' @return a integer vector. the mask, flag `i` set where `flags(i)` is `TRUE`
#' @export
bitmask_from_logical <- function(n_words, flags) {
    n_words <- .tox_as_integer_scalar(n_words, "n_words")
    flags <- .tox_as_logical_vector(flags, "flags")
    .tox_check_extent(n_words, "n_words")
    .result <- .Call("bitmask_from_logical_call", n_words, flags)
    .arguments <- c("n_bits", "n_words", "flags", "bitmask", "ierr")
    .sources <- c("flags", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$bitmask
}

#' Unpacks a bit mask into a `logical(c_bool)` array.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmask_to_logical}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags
#'   The minimum valid value is `1`.
#' @param bitmask a integer vector. the mask
#' @return a logical vector. `flags(i)` is `TRUE` where flag `i` of the mask is set
#' @export
bitmask_to_logical <- function(n_bits, bitmask) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bitmask <- .tox_as_integer_vector(bitmask, "bitmask")
    .tox_check_extent(n_bits, "n_bits")
    .result <- .Call("bitmask_to_logical_call", n_bits, bitmask)
    .arguments <- c("n_bits", "n_words", "bitmask", "flags", "ierr")
    .sources <- c(NA_character_, "bitmask", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$flags
}

#' Packs every column of a `logical(c_bool)` matrix into its own bit mask.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmasks_from_logical_2D}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_words a integer scalar. number of words per mask
#'   The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#'   The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#' @param flags a logical matrix. the flags, one mask per column
#' @return a integer matrix. the masks, column `j` packed from `flags(:, j)`
#' @export
bitmasks_from_logical_2D <- function(n_words, flags) {
    n_words <- .tox_as_integer_scalar(n_words, "n_words")
    flags <- .tox_as_logical_matrix(flags, "flags")
    .tox_check_extent(n_words, "n_words")
    .result <- .Call("bitmasks_from_logical_2D_call", n_words, flags)
    .arguments <- c("n_bits", "n_masks", "n_words", "flags", "bitmasks", "ierr")
    .sources <- c("flags", "flags", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$bitmasks
}

#' Unpacks every column of a matrix of bit masks into a `logical(c_bool)` column.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bitmasks_to_logical_2D}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags per mask
#'   The minimum valid value is `1`.
#' @param bitmasks a integer matrix. the masks, one per column
#' @return a logical matrix. the flags, column `j` unpacked from `bitmasks(:, j)`
#' @export
bitmasks_to_logical_2D <- function(n_bits, bitmasks) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bitmasks <- .tox_as_integer_matrix(bitmasks, "bitmasks")
    .tox_check_extent(n_bits, "n_bits")
    .result <- .Call("bitmasks_to_logical_2D_call", n_bits, bitmasks)
    .arguments <- c("n_bits", "n_masks", "n_words", "bitmasks", "flags", "ierr")
    .sources <- c(NA_character_, "bitmasks", "bitmasks", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$flags
}
