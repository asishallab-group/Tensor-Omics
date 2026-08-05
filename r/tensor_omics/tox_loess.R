# Generated. Do not edit.

#' Perform plain LOESS fitting
#'
#' Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
#' Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
#' response array.
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_plain}, whose argument names
#' are the ones an error message reports.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. Degree of the LOESS polynomial
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#'   The default value is `FALSE`.
#' @param save_factorization a logical scalar. Save matrix factorization flag
#'   The default value is `FALSE`.
#' @return a numeric vector. Fitted (smoothed) values of y at the evaluation points
#' @export
loess_fit_plain_expert <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence = FALSE, save_factorization = FALSE) {
    x <- .tox_as_double_vector(x, "x")
    y <- .tox_as_double_vector(y, "y")
    weights <- .tox_as_double_vector(weights, "weights")
    eval_points <- .tox_as_double_matrix(eval_points, "eval_points")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    compute_influence <- .tox_as_logical(compute_influence, "compute_influence")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    .tox_loess_required_workspace_result <- tox_loess_required_workspace(n_dim = 1L, max_neighborhood_size = max_neighborhood_size, save_factorization = save_factorization)
    int_workspace_size <- .tox_loess_required_workspace_result$int_workspace_size
    real_workspace_size <- .tox_loess_required_workspace_result$real_workspace_size

    if (length(y) != length(x))
        .tox_shape_error("y", length(y), "x", length(x))
    if (length(weights) != length(x))
        .tox_shape_error("weights", length(weights), "x", length(x))
    if (dim(eval_points)[1] != length(x))
        .tox_shape_error("eval_points", dim(eval_points)[1], "x", length(x))

    .result <- .Call("loess_fit_plain_expert_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, int_workspace_size, real_workspace_size)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "fitted_values", "ierr")
    .sources <- c("x", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "tmp_int_workspace", NA_character_, "tmp_real_workspace", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$fitted_values
}

#' Perform plain LOESS fitting
#'
#' Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
#' Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
#' response array.
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_plain_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. Degree of the LOESS polynomial
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#'   The default value is `FALSE`.
#' @param save_factorization a logical scalar. Save matrix factorization flag
#'   The default value is `FALSE`.
#' @return a numeric vector. Fitted (smoothed) values of y at the evaluation points
#' @export
loess_fit_plain <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence = FALSE, save_factorization = FALSE) {
    x <- .tox_as_double_vector(x, "x")
    y <- .tox_as_double_vector(y, "y")
    weights <- .tox_as_double_vector(weights, "weights")
    eval_points <- .tox_as_double_matrix(eval_points, "eval_points")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    compute_influence <- .tox_as_logical(compute_influence, "compute_influence")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    if (length(y) != length(x))
        .tox_shape_error("y", length(y), "x", length(x))
    if (length(weights) != length(x))
        .tox_shape_error("weights", length(weights), "x", length(x))
    if (dim(eval_points)[1] != length(x))
        .tox_shape_error("eval_points", dim(eval_points)[1], "x", length(x))

    .result <- .Call("loess_fit_plain_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "fitted_values", "ierr")
    .sources <- c("x", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$fitted_values
}

#' Perform robust LOESS fitting with bisquare reweighting
#'
#' Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
#' Fits a LOESS model to the data using robust iterations to handle outliers.
#' The robust fitting process iterates n_iters times, each iteration:
#' - Combines original weights with robust weights (down-weights from previous iteration)
#' - Runs LOESS fitting with combined weights
#' - Computes residuals (y - fitted values)
#' - Updates robust weights using bisquare function (suppresses large residuals)
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_robust}, whose argument names
#' are the ones an error message reports.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. Degree of the LOESS polynomial
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#'   The default value is `FALSE`.
#' @param save_factorization a logical scalar. Save matrix factorization flag
#'   The default value is `FALSE`.
#' @param n_iters a integer scalar. Number of robust iterations
#'   The minimum valid value is `1`.
#'   The default value is `3`.
#' @return a numeric vector. Fitted (smoothed) values of y at the evaluation points
#' @export
loess_fit_robust_expert <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence = FALSE, save_factorization = FALSE, n_iters = 3L) {
    x <- .tox_as_double_vector(x, "x")
    y <- .tox_as_double_vector(y, "y")
    weights <- .tox_as_double_vector(weights, "weights")
    eval_points <- .tox_as_double_matrix(eval_points, "eval_points")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    compute_influence <- .tox_as_logical(compute_influence, "compute_influence")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    n_iters <- .tox_as_integer_scalar(n_iters, "n_iters")
    .tox_loess_required_workspace_result <- tox_loess_required_workspace(n_dim = 1L, max_neighborhood_size = max_neighborhood_size, save_factorization = save_factorization)
    int_workspace_size <- .tox_loess_required_workspace_result$int_workspace_size
    real_workspace_size <- .tox_loess_required_workspace_result$real_workspace_size

    if (length(y) != length(x))
        .tox_shape_error("y", length(y), "x", length(x))
    if (length(weights) != length(x))
        .tox_shape_error("weights", length(weights), "x", length(x))
    if (dim(eval_points)[1] != length(x))
        .tox_shape_error("eval_points", dim(eval_points)[1], "x", length(x))

    .result <- .Call("loess_fit_robust_expert_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, n_iters, int_workspace_size, real_workspace_size)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "tmp_robust_weights", "tmp_combined_weights", "tmp_residuals", "tmp_permutation_indices", "fitted_values", "ierr")
    .sources <- c("x", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "tmp_int_workspace", NA_character_, "tmp_real_workspace", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$fitted_values
}

#' Perform robust LOESS fitting with bisquare reweighting
#'
#' Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
#' Fits a LOESS model to the data using robust iterations to handle outliers.
#' The robust fitting process iterates n_iters times, each iteration:
#' - Combines original weights with robust weights (down-weights from previous iteration)
#' - Runs LOESS fitting with combined weights
#' - Computes residuals (y - fitted values)
#' - Updates robust weights using bisquare function (suppresses large residuals)
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_robust_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. Degree of the LOESS polynomial
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#'   The default value is `FALSE`.
#' @param save_factorization a logical scalar. Save matrix factorization flag
#'   The default value is `FALSE`.
#' @param n_iters a integer scalar. Number of robust iterations
#'   The minimum valid value is `1`.
#'   The default value is `3`.
#' @return a numeric vector. Fitted (smoothed) values of y at the evaluation points
#' @export
loess_fit_robust <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence = FALSE, save_factorization = FALSE, n_iters = 3L) {
    x <- .tox_as_double_vector(x, "x")
    y <- .tox_as_double_vector(y, "y")
    weights <- .tox_as_double_vector(weights, "weights")
    eval_points <- .tox_as_double_matrix(eval_points, "eval_points")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    compute_influence <- .tox_as_logical(compute_influence, "compute_influence")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    n_iters <- .tox_as_integer_scalar(n_iters, "n_iters")
    if (length(y) != length(x))
        .tox_shape_error("y", length(y), "x", length(x))
    if (length(weights) != length(x))
        .tox_shape_error("weights", length(weights), "x", length(x))
    if (dim(eval_points)[1] != length(x))
        .tox_shape_error("eval_points", dim(eval_points)[1], "x", length(x))

    .result <- .Call("loess_fit_robust_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, n_iters)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "fitted_values", "ierr")
    .sources <- c("x", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$fitted_values
}
