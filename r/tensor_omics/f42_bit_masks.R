# Generated. Do not edit.

#' The number of words a mask of `n_bits` flags takes.
#'
#' Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_mask_n_words}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags in the mask; no flags take no words
#'   The minimum valid value is `0`.
#' @return a integer scalar. number of words the mask takes
#' @export
bit_mask_n_words <- function(n_bits) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    .result <- .Call("bit_mask_n_words_call", n_bits)
    .arguments <- c("n_bits", "n_words", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_words
}

#' Whether flag `i_bit` of a mask is set.
#'
#' Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_mask_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags in the mask
#'   The minimum valid value is `1`.
#' @param bit_mask a integer vector. the mask
#' @param i_bit a integer scalar. index of the flag, counting from 1
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_bits`.
#' @return a logical scalar. `TRUE` if the flag is set; `FALSE` for `i_bit > n_bits`, which the bits past
#'   `n_bits` would say too, as long as the mask keeps them zero
#' @export
bit_mask_test <- function(n_bits, bit_mask, i_bit) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bit_mask <- .tox_as_integer_vector(bit_mask, "bit_mask")
    i_bit <- .tox_as_integer_scalar(i_bit, "i_bit")
    .result <- .Call("bit_mask_test_call", n_bits, bit_mask, i_bit)
    .arguments <- c("n_bits", "n_words", "bit_mask", "i_bit", "is_set", "ierr")
    .sources <- c(NA_character_, "bit_mask", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_set
}

#' Packs a `logical(c_bool)` array into a bit mask.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_mask_from_logical}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_words a integer scalar. number of words of the mask
#'   The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#'   The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#' @param flags a logical vector. the flags to pack
#' @return a integer vector. the mask, flag `i` set where `flags(i)` is `TRUE`
#' @export
bit_mask_from_logical <- function(n_words, flags) {
    n_words <- .tox_as_integer_scalar(n_words, "n_words")
    flags <- .tox_as_logical_vector(flags, "flags")
    .result <- .Call("bit_mask_from_logical_call", n_words, flags)
    .arguments <- c("n_bits", "n_words", "flags", "bit_mask", "ierr")
    .sources <- c("flags", "bit_mask", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$bit_mask
}

#' Unpacks a bit mask into a `logical(c_bool)` array.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_mask_to_logical}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags
#'   The minimum valid value is `1`.
#' @param bit_mask a integer vector. the mask
#' @return a logical vector. `flags(i)` is `TRUE` where flag `i` of the mask is set
#' @export
bit_mask_to_logical <- function(n_bits, bit_mask) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bit_mask <- .tox_as_integer_vector(bit_mask, "bit_mask")
    .result <- .Call("bit_mask_to_logical_call", n_bits, bit_mask)
    .arguments <- c("n_bits", "n_words", "bit_mask", "flags", "ierr")
    .sources <- c("flags", "bit_mask", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$flags
}

#' Packs every column of a `logical(c_bool)` matrix into its own bit mask.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_masks_from_logical_2D}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_words a integer scalar. number of words per mask
#'   The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#'   The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
#' @param flags a logical matrix. the flags, one mask per column
#' @return a integer matrix. the masks, column `j` packed from `flags(:, j)`
#' @export
bit_masks_from_logical_2D <- function(n_words, flags) {
    n_words <- .tox_as_integer_scalar(n_words, "n_words")
    flags <- .tox_as_logical_matrix(flags, "flags")
    .result <- .Call("bit_masks_from_logical_2D_call", n_words, flags)
    .arguments <- c("n_bits", "n_masks", "n_words", "flags", "bit_masks", "ierr")
    .sources <- c("flags", "flags", "bit_masks", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$bit_masks
}

#' Unpacks every column of a matrix of bit masks into a `logical(c_bool)` column.
#'
#' Generated from the Fortran procedure \code{f42_bit_masks::bit_masks_to_logical_2D}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bits a integer scalar. number of flags per mask
#'   The minimum valid value is `1`.
#' @param bit_masks a integer matrix. the masks, one per column
#' @return a logical matrix. the flags, column `j` unpacked from `bit_masks(:, j)`
#' @export
bit_masks_to_logical_2D <- function(n_bits, bit_masks) {
    n_bits <- .tox_as_integer_scalar(n_bits, "n_bits")
    bit_masks <- .tox_as_integer_matrix(bit_masks, "bit_masks")
    .result <- .Call("bit_masks_to_logical_2D_call", n_bits, bit_masks)
    .arguments <- c("n_bits", "n_masks", "n_words", "bit_masks", "flags", "ierr")
    .sources <- c("flags", "bit_masks", "bit_masks", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$flags
}
