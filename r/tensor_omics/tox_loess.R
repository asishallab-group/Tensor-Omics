# Generated. Do not edit.

#' Recommend workspace sizes based on Netlib exact formulas
#'
#' Computes the required sizes for integer and real workspace arrays.
#' These sizes depend on the dimensionality of the data and the maximum neighborhood size.
#'
#' @param n_dim a integer scalar. Dimensionality of the data
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param save_factorization a logical scalar. Save matrix factorization flag
#' @return a named list with elements `int_workspace_size`, `real_workspace_size`.
#'
#' Generated from the Fortran procedure \code{tox_loess::tox_loess_required_workspace}.
#' @export
tox_loess_required_workspace <- function(n_dim, max_neighborhood_size, save_factorization) {
    n_dim <- .tox_as_integer_scalar(n_dim, "n_dim")
    max_neighborhood_size <- .tox_as_integer_scalar(max_neighborhood_size, "max_neighborhood_size")
    save_factorization <- .tox_as_logical(save_factorization, "save_factorization")
    .result <- .Call("tox_loess_required_workspace_call", n_dim, max_neighborhood_size, save_factorization)
    .arguments <- c("n_dim", "max_neighborhood_size", "int_workspace_size", "real_workspace_size", "save_factorization")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        int_workspace_size = .result$int_workspace_size,
        real_workspace_size = .result$real_workspace_size
    )
}

#' Perform plain LOESS fitting
#'
#' Fits a LOESS model to the data using the specified smoothing parameter.
#' Outputs the smoothed response variable array. Caller-provided workspace must already be
#' sized via [[tox_loess(module):tox_loess_required_workspace(subroutine)]]; no input validation
#' beyond `n`, `span`, `degree`, and workspace-size bounds is performed here.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#' @param degree a integer scalar. Degree of the LOESS polynomial
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#' @param save_factorization a logical scalar. Save matrix factorization flag
#' @return Fitted (smoothed) values of y at the evaluation points
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_plain}.
#' @export
loess_fit_plain <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization) {
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

    .result <- .Call("loess_fit_plain_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, int_workspace_size, real_workspace_size)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "fitted_values", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$fitted_values
}

#' Perform robust LOESS fitting with bisquare reweighting
#'
#' Fits a LOESS model to the data using robust iterations to handle outliers.
#' The robust fitting process iterates n_iters times, each iteration:
#' - Combines original weights with robust weights (down-weights from previous iteration)
#' - Runs LOESS fitting with combined weights
#' - Computes residuals (y - fitted values)
#' - Updates robust weights using bisquare function (suppresses large residuals)
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param weights a numeric vector. Weight array for data points
#' @param eval_points a numeric matrix. Evaluation points (x values at which the fitted curve is computed)
#' @param span a numeric scalar. Smoothing parameter for LOESS
#' @param degree a integer scalar. Degree of the LOESS polynomial
#' @param max_neighborhood_size a integer scalar. Maximum neighborhood size
#' @param compute_influence a logical scalar. Influence calculation flag
#' @param save_factorization a logical scalar. Save matrix factorization flag
#' @param n_iters a integer scalar. Number of robust iterations
#' @param robust_weights a numeric vector. Robust bisquare weights (updated each iteration, initialized to 1.0)
#' @param combined_weights a numeric vector. Combined weights: product of user weights and robust weights (weights(i) * robust_weights(i))
#' @param residuals a numeric vector. Residuals (y - fitted_values), used to compute bisquare robust weights
#' @param permutation_indices a integer vector. Permutation indices array (from NetLib bisquare weight computation)
#' @return a named list with elements `robust_weights`, `combined_weights`, `residuals`, `permutation_indices`, `fitted_values`.
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_fit_robust}.
#' @export
loess_fit_robust <- function(x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, n_iters, robust_weights, combined_weights, residuals, permutation_indices) {
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
    robust_weights <- .tox_as_double_vector(robust_weights, "robust_weights")
    combined_weights <- .tox_as_double_vector(combined_weights, "combined_weights")
    residuals <- .tox_as_double_vector(residuals, "residuals")
    permutation_indices <- .tox_as_integer_vector(permutation_indices, "permutation_indices")
    .tox_loess_required_workspace_result <- tox_loess_required_workspace(n_dim = 1L, max_neighborhood_size = max_neighborhood_size, save_factorization = save_factorization)
    int_workspace_size <- .tox_loess_required_workspace_result$int_workspace_size
    real_workspace_size <- .tox_loess_required_workspace_result$real_workspace_size

    if (length(y) != length(x))
        .tox_shape_error("y", length(y), "x", length(x))
    if (length(weights) != length(x))
        .tox_shape_error("weights", length(weights), "x", length(x))
    if (dim(eval_points)[1] != length(x))
        .tox_shape_error("eval_points", dim(eval_points)[1], "x", length(x))
    if (length(robust_weights) != length(x))
        .tox_shape_error("robust_weights", length(robust_weights), "x", length(x))
    if (length(combined_weights) != length(x))
        .tox_shape_error("combined_weights", length(combined_weights), "x", length(x))
    if (length(residuals) != length(x))
        .tox_shape_error("residuals", length(residuals), "x", length(x))
    if (length(permutation_indices) != length(x))
        .tox_shape_error("permutation_indices", length(permutation_indices), "x", length(x))

    .result <- .Call("loess_fit_robust_call", x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, n_iters, int_workspace_size, real_workspace_size, robust_weights, combined_weights, residuals, permutation_indices)
    .arguments <- c("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "robust_weights", "combined_weights", "residuals", "permutation_indices", "fitted_values", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        robust_weights = .result$robust_weights,
        combined_weights = .result$combined_weights,
        residuals = .result$residuals,
        permutation_indices = .result$permutation_indices,
        fitted_values = .result$fitted_values
    )
}

#' Wrapper subroutine for LOESS fitting (plain or robust)
#'
#' This subroutine selects between plain and robust LOESS fitting based on the mode.
#' It dynamically allocates the required arrays and computes workspace sizes, and handles
#' degenerate inputs (single point, near-constant `x`, or fewer unique `x` values than the
#' polynomial degree requires) by falling back to an identity/copy mapping instead of calling
#' into netlib.
#'
#' Parameters:
#' - mode: Specifies the type of LOESS fitting to perform.
#' - 0: Plain LOESS fitting. This mode performs a single pass of LOESS fitting without any additional weighting or iterations. It is suitable for datasets without significant outliers.
#' - 1: Robust LOESS fitting. This mode applies bisquare reweighting over multiple iterations to reduce the influence of outliers. The number of iterations is controlled by the `n_iters` parameter.
#'
#' @param x a numeric vector. Predictor variable array
#' @param y a numeric vector. Response variable array
#' @param span a numeric scalar. Smoothing parameter for LOESS
#' @param degree a integer scalar. Degree of the LOESS polynomial
#' @param mode a string, one of "plain", "robust". Mode of operation
#' @param n_iters a integer scalar. Number of robust iterations, ignored in [[tox_loess(module):MODE_PLAIN(variable)]].
#' @return Fitted (smoothed) values of y
#'
#' Generated from the Fortran procedure \code{tox_loess::loess_alloc}.
#' @export
loess <- function(x, y, span, degree, mode, n_iters = 3L) {
    x <- .tox_as_double_vector(x, "x")
    y <- .tox_as_double_vector(y, "y")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    mode <- .tox_as_mode(mode, "mode", c("plain", "robust"))
    n_iters <- .tox_as_integer_scalar(n_iters, "n_iters")
    .result <- .Call("loess_call", x, y, span, degree, mode, n_iters)
    .arguments <- c("x", "y", "span", "degree", "fitted_values", "mode", "n_iters", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$fitted_values
}
