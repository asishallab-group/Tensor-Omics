#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_loess(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_loess_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: loess_fit_plain_expert_c
    public :: loess_fit_plain_c
    public :: loess_fit_robust_expert_c
    public :: loess_fit_robust_c
    public :: loess_c

contains

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_plain(subroutine)]]
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array. Caller-provided workspace must already be sized via
    !| [[tox_loess_kernel(module):tox_loess_required_workspace(subroutine)]].
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
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
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
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
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
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array. Caller-provided workspace must already be sized via
    !| [[tox_loess_kernel(module):tox_loess_required_workspace(subroutine)]].
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
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
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
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
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
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
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
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
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
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n), intent(in), target :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n, 1), intent(in), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
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
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
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

    !> summary: C-wrapper for [[tox_loess(module):loess_alloc(subroutine)]]
    !| This kernel selects between plain and robust LOESS fitting based on the mode. It dynamically
    !| allocates the required arrays and computes workspace sizes, and handles degenerate inputs (single
    !| point, near-constant `x`, or fewer unique `x` values than the polynomial degree requires) by
    !| falling back to an identity/copy mapping instead of calling into netlib. As a self-allocating
    !| pipeline it validates its own inputs, so the generated wrapper adds only the binding surface.
    !|
    !| Parameters:
    !| - mode: Specifies the type of LOESS fitting to perform.
    !| - 0: Plain LOESS fitting. This mode performs a single pass of LOESS fitting without any additional weighting or iterations. It is suitable for datasets without significant outliers.
    !| - 1: Robust LOESS fitting. This mode applies bisquare reweighting over multiple iterations to reduce the influence of outliers. The number of iterations is controlled by the `n_iters` parameter.
    subroutine loess_c(&
            x,&
            n_x_elements,&
            y,&
            n_y_elements,&
            span,&
            degree,&
            fitted_values,&
            mode,&
            n_iters,&
            ierr&
        ) bind(C, name="loess_c")
        use tox_loess, only: loess_alloc
        use tox_loess_kernel, only: MODE_PLAIN, MODE_ROBUST

        integer(c_int), intent(in), target :: n_x_elements
            !! number of elements in `x`
        real(c_double), dimension(n_y_elements), intent(in), target :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), intent(in), target :: n_y_elements
            !! number of elements in `y`
        real(c_double), dimension(n_x_elements), intent(in), target :: x
            !! Predictor variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
        real(c_double), dimension(size(y)), intent(out), target :: fitted_values
            !! Fitted (smoothed) values of y
        character(len=1, kind=c_char), dimension(6), intent(in), target :: mode
            !! Mode of operation
            !!
            !! | Mode                 | Value                                              |
            !! |----------------------|----------------------------------------------------|
            !! | Plain LOESS fitting  | [[tox_loess_kernel(module):MODE_PLAIN(variable)]]  |
            !! | Robust LOESS fitting | [[tox_loess_kernel(module):MODE_ROBUST(variable)]] |
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations, ignored in [[tox_loess_kernel(module):MODE_PLAIN(variable)]].
            !! The default value is `3_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_x_elements)
        M_CHECK_NON_NULL(n_y_elements)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_ARRAY_NON_NULL(x, n_x_elements)
        M_CHECK_ARRAY_NON_NULL(y, n_y_elements)
        M_CHECK_ARRAY_NON_NULL(fitted_values, (size(y)))
        M_CHECK_ARRAY_NON_NULL(mode, 6)

        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("plain")
                    mode_mode_f = MODE_PLAIN
                case ("robust")
                    mode_mode_f = MODE_ROBUST
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call loess_alloc(&
            x = x,&
            y = y,&
            span = span,&
            degree = degree,&
            fitted_values = fitted_values,&
            mode = mode_mode_f,&
            n_iters = n_iters,&
            ierr = ierr&
        )
    end subroutine loess_c

end module tox_loess_c
#endif
