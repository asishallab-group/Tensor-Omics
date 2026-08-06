#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_loess(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_loess_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: loess_fit_plain_expert_c
    public :: loess_fit_plain_c
    public :: loess_fit_robust_expert_c
    public :: loess_fit_robust_c

contains

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_plain(subroutine)]]
    !| Data too degenerate to fit is answered directly, by the observations themselves; see
    !| [[tox_loess_kernel(module):loess_degenerate_fit]].
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array.
    subroutine loess_fit_plain_expert_c(&
            n,&
            x,&
            y,&
            weights,&
            eval_points,&
            span,&
            degree,&
            max_neighborhood_size,&
            compute_influence,&
            save_factorization,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_hat_diag,&
            fitted_values,&
            ierr&
        ) bind(C, name="loess_fit_plain_expert_c")
        use tox_loess, only: loess_fit_plain

        integer(c_int), intent(in), target :: n
            !! Total number of data points
        integer(c_int), intent(in), target :: int_workspace_size
            !! Required size of the integer workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `10000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Required size of the real workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `100000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        real(c_double), dimension(n), intent(in), target :: x
            !! Predictor variable array
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(c_int), intent(in), target :: max_neighborhood_size
            !! Maximum neighborhood size
        logical(c_bool), intent(in), target :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n), intent(out), target :: tmp_hat_diag
            !! Diagonal elements of the hat matrix
        real(c_double), dimension(n), intent(out), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: compute_influence_f
        logical :: save_factorization_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(max_neighborhood_size)
        M_CHECK_NON_NULL(compute_influence)
        M_CHECK_NON_NULL(save_factorization)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(x, n)
        M_CHECK_ARRAY_NON_NULL(y, n)
        M_CHECK_ARRAY_NON_NULL(weights, n)
        M_CHECK_ARRAY_NON_NULL(eval_points, n * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_hat_diag, n)
        M_CHECK_ARRAY_NON_NULL(fitted_values, n)

        compute_influence_f = compute_influence
        save_factorization_f = save_factorization

        call loess_fit_plain(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence_f,&
            save_factorization = save_factorization_f,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_plain_expert_c

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_plain_alloc(subroutine)]]
    !| Data too degenerate to fit is answered directly, by the observations themselves; see
    !| [[tox_loess_kernel(module):loess_degenerate_fit]].
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array.
    subroutine loess_fit_plain_c(&
            n,&
            x,&
            y,&
            weights,&
            eval_points,&
            span,&
            degree,&
            max_neighborhood_size,&
            compute_influence,&
            save_factorization,&
            fitted_values,&
            ierr&
        ) bind(C, name="loess_fit_plain_c")
        use tox_loess, only: loess_fit_plain_alloc

        integer(c_int), intent(in), target :: n
            !! Total number of data points
        real(c_double), dimension(n), intent(in), target :: x
            !! Predictor variable array
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(c_int), intent(in), target :: max_neighborhood_size
            !! Maximum neighborhood size
        logical(c_bool), intent(in), target :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        real(c_double), dimension(n), intent(out), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: compute_influence_f
        logical :: save_factorization_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(max_neighborhood_size)
        M_CHECK_NON_NULL(compute_influence)
        M_CHECK_NON_NULL(save_factorization)
        M_CHECK_ARRAY_NON_NULL(x, n)
        M_CHECK_ARRAY_NON_NULL(y, n)
        M_CHECK_ARRAY_NON_NULL(weights, n)
        M_CHECK_ARRAY_NON_NULL(eval_points, n * 1)
        M_CHECK_ARRAY_NON_NULL(fitted_values, n)

        compute_influence_f = compute_influence
        save_factorization_f = save_factorization

        call loess_fit_plain_alloc(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence_f,&
            save_factorization = save_factorization_f,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_plain_c

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_robust(subroutine)]]
    !| Data too degenerate to fit is answered directly, by the observations themselves; see
    !| [[tox_loess_kernel(module):loess_degenerate_fit]].
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !| - Combines original weights with robust weights (down-weights from previous iteration)
    !| - Runs LOESS fitting with combined weights
    !| - Computes residuals (y - fitted values)
    !| - Updates robust weights using bisquare function (suppresses large residuals)
    subroutine loess_fit_robust_expert_c(&
            n,&
            x,&
            y,&
            weights,&
            eval_points,&
            span,&
            degree,&
            max_neighborhood_size,&
            compute_influence,&
            save_factorization,&
            n_iters,&
            tmp_int_workspace,&
            int_workspace_size,&
            tmp_real_workspace,&
            real_workspace_size,&
            tmp_hat_diag,&
            tmp_robust_weights,&
            tmp_combined_weights,&
            tmp_residuals,&
            tmp_permutation_indices,&
            fitted_values,&
            ierr&
        ) bind(C, name="loess_fit_robust_expert_c")
        use tox_loess, only: loess_fit_robust

        integer(c_int), intent(in), target :: n
            !! Total number of data points
        integer(c_int), intent(in), target :: int_workspace_size
            !! Required size of the integer workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `10000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        integer(c_int), intent(in), target :: real_workspace_size
            !! Required size of the real workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `100000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        real(c_double), dimension(n), intent(in), target :: x
            !! Predictor variable array
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(c_int), intent(in), target :: max_neighborhood_size
            !! Maximum neighborhood size
        logical(c_bool), intent(in), target :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
            !! The default value is `3_int32`.
        integer(c_int), dimension(int_workspace_size), intent(out), target :: tmp_int_workspace
            !! Integer workspace array
        real(c_double), dimension(real_workspace_size), intent(out), target :: tmp_real_workspace
            !! Real workspace array
        real(c_double), dimension(n), intent(out), target :: tmp_hat_diag
            !! Diagonal elements of the hat matrix
        real(c_double), dimension(n), intent(out), target :: tmp_robust_weights
            !! Robust bisquare weights (updated each iteration, initialized to 1.0)
        real(c_double), dimension(n), intent(out), target :: tmp_combined_weights
            !! Combined weights: product of user weights and robust weights (weights(i) * robust_weights(i))
        real(c_double), dimension(n), intent(out), target :: tmp_residuals
            !! Residuals (y - fitted_values), used to compute bisquare robust weights
        integer(c_int), dimension(n), intent(out), target :: tmp_permutation_indices
            !! Permutation indices array (from NetLib bisquare weight computation)
        real(c_double), dimension(n), intent(out), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: compute_influence_f
        logical :: save_factorization_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(max_neighborhood_size)
        M_CHECK_NON_NULL(compute_influence)
        M_CHECK_NON_NULL(save_factorization)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(x, n)
        M_CHECK_ARRAY_NON_NULL(y, n)
        M_CHECK_ARRAY_NON_NULL(weights, n)
        M_CHECK_ARRAY_NON_NULL(eval_points, n * 1)
        M_CHECK_ARRAY_NON_NULL(tmp_int_workspace, int_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_real_workspace, real_workspace_size)
        M_CHECK_ARRAY_NON_NULL(tmp_hat_diag, n)
        M_CHECK_ARRAY_NON_NULL(tmp_robust_weights, n)
        M_CHECK_ARRAY_NON_NULL(tmp_combined_weights, n)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals, n)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_indices, n)
        M_CHECK_ARRAY_NON_NULL(fitted_values, n)

        compute_influence_f = compute_influence
        save_factorization_f = save_factorization

        call loess_fit_robust(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence_f,&
            save_factorization = save_factorization_f,&
            n_iters = n_iters,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            tmp_robust_weights = tmp_robust_weights,&
            tmp_combined_weights = tmp_combined_weights,&
            tmp_residuals = tmp_residuals,&
            tmp_permutation_indices = tmp_permutation_indices,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_robust_expert_c

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_robust_alloc(subroutine)]]
    !| Data too degenerate to fit is answered directly, by the observations themselves; see
    !| [[tox_loess_kernel(module):loess_degenerate_fit]].
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !| - Combines original weights with robust weights (down-weights from previous iteration)
    !| - Runs LOESS fitting with combined weights
    !| - Computes residuals (y - fitted values)
    !| - Updates robust weights using bisquare function (suppresses large residuals)
    subroutine loess_fit_robust_c(&
            n,&
            x,&
            y,&
            weights,&
            eval_points,&
            span,&
            degree,&
            max_neighborhood_size,&
            compute_influence,&
            save_factorization,&
            n_iters,&
            fitted_values,&
            ierr&
        ) bind(C, name="loess_fit_robust_c")
        use tox_loess, only: loess_fit_robust_alloc

        integer(c_int), intent(in), target :: n
            !! Total number of data points
        real(c_double), dimension(n), intent(in), target :: x
            !! Predictor variable array
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(c_int), intent(in), target :: max_neighborhood_size
            !! Maximum neighborhood size
        logical(c_bool), intent(in), target :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
            !! The default value is `3_int32`.
        real(c_double), dimension(n), intent(out), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: compute_influence_f
        logical :: save_factorization_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(max_neighborhood_size)
        M_CHECK_NON_NULL(compute_influence)
        M_CHECK_NON_NULL(save_factorization)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_ARRAY_NON_NULL(x, n)
        M_CHECK_ARRAY_NON_NULL(y, n)
        M_CHECK_ARRAY_NON_NULL(weights, n)
        M_CHECK_ARRAY_NON_NULL(eval_points, n * 1)
        M_CHECK_ARRAY_NON_NULL(fitted_values, n)

        compute_influence_f = compute_influence
        save_factorization_f = save_factorization

        call loess_fit_robust_alloc(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence_f,&
            save_factorization = save_factorization_f,&
            n_iters = n_iters,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_robust_c

end module tox_loess_c
#endif
