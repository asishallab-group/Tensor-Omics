# Generated. Do not edit.

#' Whether a grown ensemble at t+1 is still compatible with its own growth trajectory
#'
#' Four criteria, all must hold: (1) tangent-space drift -- for every reference basis in
#' {U_first} union {U_history(:,:,1:history_len)} that shares d_tp1's rank, the chordal
#' distance to U_tp1 (via `dgesvd` on M = U_r(:,1:d)^T U_tp1(:,1:d), whose singular values
#' are cos(theta_i) directly) must not exceed
#' chordal_dist_max_as_prcnt_of_range * sqrt(d_tp1); references with a different d are
#' skipped for this criterion (no shared dimension to compare angles over), judged instead
#' by criterion (2). Only the candidate is compared against each reference, not every pair
#' within the reference set itself -- sufficient, not a shortcut, since every reference
#' already passed this same criterion against every other reference present in the set at
#' the growth step it was itself accepted (an inductive invariant FIFO eviction cannot
#' break), so this costs O(history_len) small SVDs per growth step, not O(history_len^2).
#' (2) max(|d_tp1 - d_first|, |d_tp1 - d_history(history_len)|) <= d_max. (3)
#' |log(G_tp1/G_t)| <= G_max, where G_t = the most recently accepted iteration's spectral
#' gap. (4) |log(RMSE_tp1/RMSE_t)| <= RMSE_change_max, where RMSE = sqrt(normal_error) --
#' see `normal_error` in the sibling observable kernel module; already free, no new SVD or
#' storage. `ierr` is set only if a LAPACK SVD fails to converge -- not a condition any
#' input check could foresee.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_accept::accept_ensemble}, whose argument names
#' are the ones an error message reports.
#'
#' @param U_first a numeric matrix. Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
#'   `misc/mod_STC.md`, "Output"
#' @param d_first a integer scalar. Ensemble's intrinsic dimension at its bootstrap iteration
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param U_history a numeric array of rank 3. Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
#'   valid (see `history_len`)
#' @param d_history a integer vector. Trailing intrinsic dimensions, one per column of U_history; only entries
#'   1:history_len are valid
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param history_len a integer scalar. Number of valid columns in U_history/d_history; column history_len is the most
#'   recently accepted iteration
#'   The minimum valid value is `1`.
#'   The maximum valid value is `o`.
#' @param G_t a numeric scalar. Ensemble's spectral gap at the most recently accepted iteration
#'   The minimum valid value is `above(0.0)`.
#' @param normal_error_t a numeric scalar. Ensemble's normal_error at the most recently accepted iteration, see
#'   `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
#'   noise at all in its normal directions -- unlike G_t, which is a ratio already
#'   protected by its own +epsilon denominator (see `observable`) and so is
#'   required to be strictly positive.
#'   The minimum valid value is `0.0`.
#' @param U_tp1 a numeric matrix. Candidate ensemble's tangent+normal basis
#' @param d_tp1 a integer scalar. Candidate ensemble's intrinsic dimension
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_dimensions`.
#' @param G_tp1 a numeric scalar. Candidate ensemble's spectral gap
#'   The minimum valid value is `above(0.0)`.
#' @param normal_error_tp1 a numeric scalar. Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
#'   The minimum valid value is `0.0`.
#' @param chordal_dist_max_as_prcnt_of_range a numeric scalar. Maximum tolerated chordal distance between tangent bases, as a fraction of its
#'   own [0, sqrt(d)] range, see `accept_ensemble`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param d_max a integer scalar. Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
#'   The minimum valid value is `0`.
#' @param G_max a numeric scalar. Maximum tolerated |log(G_tp1/G_t)|
#'   The minimum valid value is `0.0`.
#' @param RMSE_change_max a numeric scalar. Maximum tolerated |log(RMSE_tp1/RMSE_t)|
#'   The minimum valid value is `0.0`.
#' @return a logical scalar. TRUE if all four acceptance criteria are satisfied
#' @export
accept_ensemble <- function(U_first, d_first, U_history, d_history, history_len, G_t, normal_error_t, U_tp1, d_tp1, G_tp1, normal_error_tp1, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max) {
    U_first <- .tox_as_double_matrix(U_first, "U_first")
    d_first <- .tox_as_integer_scalar(d_first, "d_first")
    U_history <- .tox_as_double_array(U_history, "U_history", 3L)
    d_history <- .tox_as_integer_vector(d_history, "d_history")
    history_len <- .tox_as_integer_scalar(history_len, "history_len")
    G_t <- .tox_as_double_scalar(G_t, "G_t")
    normal_error_t <- .tox_as_double_scalar(normal_error_t, "normal_error_t")
    U_tp1 <- .tox_as_double_matrix(U_tp1, "U_tp1")
    d_tp1 <- .tox_as_integer_scalar(d_tp1, "d_tp1")
    G_tp1 <- .tox_as_double_scalar(G_tp1, "G_tp1")
    normal_error_tp1 <- .tox_as_double_scalar(normal_error_tp1, "normal_error_tp1")
    chordal_dist_max_as_prcnt_of_range <- .tox_as_double_scalar(chordal_dist_max_as_prcnt_of_range, "chordal_dist_max_as_prcnt_of_range")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    RMSE_change_max <- .tox_as_double_scalar(RMSE_change_max, "RMSE_change_max")
    if (dim(U_history)[1] != dim(U_first)[1])
        .tox_shape_error("U_history", dim(U_history)[1], "U_first", dim(U_first)[1])
    if (dim(U_tp1)[1] != dim(U_first)[1])
        .tox_shape_error("U_tp1", dim(U_tp1)[1], "U_first", dim(U_first)[1])
    if (length(d_history) != dim(U_history)[3])
        .tox_shape_error("d_history", length(d_history), "U_history", dim(U_history)[3])

    .result <- .Call("accept_ensemble_call", U_first, d_first, U_history, d_history, history_len, G_t, normal_error_t, U_tp1, d_tp1, G_tp1, normal_error_tp1, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max)
    .arguments <- c("n_dimensions", "o", "U_first", "d_first", "U_history", "d_history", "history_len", "G_t", "normal_error_t", "U_tp1", "d_tp1", "G_tp1", "normal_error_tp1", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "is_accepted", "ierr")
    .sources <- c("U_first", "U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_accepted
}
