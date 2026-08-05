#include <src/macros.h>

!> summary: Wrappers for [[tox_loess_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_loess
    use tox_loess_kernel, only: EPS_LOESS, loess_degenerate_fit, loess_fit_plain_kernel, loess_fit_robust_kernel
    use tox_loess_kernel, only: tox_loess_required_workspace
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: loess_fit_plain
    public :: loess_fit_plain_alloc
    public :: loess_fit_robust
    public :: loess_fit_robust_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_loess_kernel(module):loess_fit_plain_kernel]].
    !| Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array. Caller-provided workspace must already be sized via
    !| [[tox_loess_kernel(module):tox_loess_required_workspace(subroutine)]].
    subroutine loess_fit_plain(&
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
        )
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `10000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `100000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        real(real64), dimension(n), intent(in) :: x
            !! Predictor variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n, 1), intent(in) :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n), intent(out) :: tmp_hat_diag
            !! Diagonal elements of the hat matrix
        real(real64), dimension(n), intent(out) :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        logical :: handled

        call set_ok(ierr)
        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(span, ierr, arg_pos=6_int32, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(degree, ierr, arg_pos=7_int32, min=0_int32, max=2_int32)
        call validate_in_range_int(int_workspace_size, ierr, arg_pos=12_int32, min=10000_int32)
        call validate_in_range_int(real_workspace_size, ierr, arg_pos=14_int32, min=100000_int32)
        if (is_err(ierr)) return

        call loess_degenerate_fit(&
            n = n,&
            x = x,&
            y = y,&
            degree = degree,&
            fitted_values = fitted_values,&
            handled = handled,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        if (handled) return

        call loess_fit_plain_kernel(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence,&
            save_factorization = save_factorization,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_plain

    !> summary: Allocates its work arrays, then calls [[tox_loess_kernel(module):loess_fit_plain_kernel]].
    !| Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array. Caller-provided workspace must already be sized via
    !| [[tox_loess_kernel(module):tox_loess_required_workspace(subroutine)]].
    subroutine loess_fit_plain_alloc(&
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
        )
        integer(int32), intent(in) :: n
            !! Total number of data points
        real(real64), dimension(n), intent(in) :: x
            !! Predictor variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n, 1), intent(in) :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        real(real64), dimension(n), intent(out) :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_hat_diag
        logical :: handled

        call set_ok(ierr)
        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(span, ierr, arg_pos=6_int32, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(degree, ierr, arg_pos=7_int32, min=0_int32, max=2_int32)
        if (is_err(ierr)) return

        call loess_degenerate_fit(&
            n = n,&
            x = x,&
            y = y,&
            degree = degree,&
            fitted_values = fitted_values,&
            handled = handled,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        if (handled) return

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = max_neighborhood_size,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = save_factorization&
        )
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_hat_diag(n))

        call loess_fit_plain_kernel(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence,&
            save_factorization = save_factorization,&
            tmp_int_workspace = tmp_int_workspace,&
            int_workspace_size = int_workspace_size,&
            tmp_real_workspace = tmp_real_workspace,&
            real_workspace_size = real_workspace_size,&
            tmp_hat_diag = tmp_hat_diag,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )
    end subroutine loess_fit_plain_alloc

    !> summary: Validates its inputs, then calls [[tox_loess_kernel(module):loess_fit_robust_kernel]].
    !| Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !| - Combines original weights with robust weights (down-weights from previous iteration)
    !| - Runs LOESS fitting with combined weights
    !| - Computes residuals (y - fitted values)
    !| - Updates robust weights using bisquare function (suppresses large residuals)
    subroutine loess_fit_robust(&
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
        )
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `int_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `10000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array
            !! It is *VERY IMPORTANT* to compute this argument from the `real_workspace_size` output produced by [[tox_loess_kernel(module):tox_loess_required_workspace]].
            !! The minimum valid value is `100000_int32`.
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_dim          | 1_int32     |
        real(real64), dimension(n), intent(in) :: x
            !! Predictor variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n, 1), intent(in) :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
            !! The default value is `3_int32`.
        integer(int32), dimension(int_workspace_size), intent(out) :: tmp_int_workspace
            !! Integer workspace array
        real(real64), dimension(real_workspace_size), intent(out) :: tmp_real_workspace
            !! Real workspace array
        real(real64), dimension(n), intent(out) :: tmp_hat_diag
            !! Diagonal elements of the hat matrix
        real(real64), dimension(n), intent(out) :: tmp_robust_weights
            !! Robust bisquare weights (updated each iteration, initialized to 1.0)
        real(real64), dimension(n), intent(out) :: tmp_combined_weights
            !! Combined weights: product of user weights and robust weights (weights(i) * robust_weights(i))
        real(real64), dimension(n), intent(out) :: tmp_residuals
            !! Residuals (y - fitted_values), used to compute bisquare robust weights
        integer(int32), dimension(n), intent(out) :: tmp_permutation_indices
            !! Permutation indices array (from NetLib bisquare weight computation)
        real(real64), dimension(n), intent(out) :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        logical :: handled

        call set_ok(ierr)
        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(span, ierr, arg_pos=6_int32, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(degree, ierr, arg_pos=7_int32, min=0_int32, max=2_int32)
        call validate_in_range_int(n_iters, ierr, arg_pos=11_int32, min=1_int32)
        call validate_in_range_int(int_workspace_size, ierr, arg_pos=13_int32, min=10000_int32)
        call validate_in_range_int(real_workspace_size, ierr, arg_pos=15_int32, min=100000_int32)
        if (is_err(ierr)) return

        call loess_degenerate_fit(&
            n = n,&
            x = x,&
            y = y,&
            degree = degree,&
            fitted_values = fitted_values,&
            handled = handled,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        if (handled) return

        call loess_fit_robust_kernel(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence,&
            save_factorization = save_factorization,&
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
    end subroutine loess_fit_robust

    !> summary: Allocates its work arrays, then calls [[tox_loess_kernel(module):loess_fit_robust_kernel]].
    !| Every generated wrapper runs [[tox_loess_kernel(module):loess_degenerate_fit]] first, which may handle the call and skip this one.
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !| - Combines original weights with robust weights (down-weights from previous iteration)
    !| - Runs LOESS fitting with combined weights
    !| - Computes residuals (y - fitted values)
    !| - Updates robust weights using bisquare function (suppresses large residuals)
    subroutine loess_fit_robust_alloc(&
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
        )
        integer(int32), intent(in) :: n
            !! Total number of data points
        real(real64), dimension(n), intent(in) :: x
            !! Predictor variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: y
            !! Response variable array
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n), intent(in) :: weights
            !! Weight array for data points
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n, 1), intent(in) :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! The minimum valid value is `EPS_LOESS`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `2_int32`.
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! The default value is `.false.`.
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! The default value is `.false.`.
        integer(int32), intent(in), optional :: n_iters
            !! Number of robust iterations
            !! The minimum valid value is `1_int32`.
            !! The default value is `3_int32`.
        real(real64), dimension(n), intent(out) :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_int_workspace
        integer(int32) :: int_workspace_size
        real(real64), dimension(:), allocatable :: tmp_real_workspace
        integer(int32) :: real_workspace_size
        real(real64), dimension(:), allocatable :: tmp_hat_diag
        real(real64), dimension(:), allocatable :: tmp_robust_weights
        real(real64), dimension(:), allocatable :: tmp_combined_weights
        real(real64), dimension(:), allocatable :: tmp_residuals
        integer(int32), dimension(:), allocatable :: tmp_permutation_indices
        logical :: handled

        call set_ok(ierr)
        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(span, ierr, arg_pos=6_int32, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(degree, ierr, arg_pos=7_int32, min=0_int32, max=2_int32)
        call validate_in_range_int(n_iters, ierr, arg_pos=11_int32, min=1_int32)
        if (is_err(ierr)) return

        call loess_degenerate_fit(&
            n = n,&
            x = x,&
            y = y,&
            degree = degree,&
            fitted_values = fitted_values,&
            handled = handled,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        if (handled) return

        call tox_loess_required_workspace(&
            n_dim = 1_int32,&
            max_neighborhood_size = max_neighborhood_size,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = save_factorization&
        )
        M_ALLOCATE(tmp_int_workspace(int_workspace_size))
        M_ALLOCATE(tmp_real_workspace(real_workspace_size))
        M_ALLOCATE(tmp_hat_diag(n))
        M_ALLOCATE(tmp_robust_weights(n))
        M_ALLOCATE(tmp_combined_weights(n))
        M_ALLOCATE(tmp_residuals(n))
        M_ALLOCATE(tmp_permutation_indices(n))

        call loess_fit_robust_kernel(&
            n = n,&
            x = x,&
            y = y,&
            weights = weights,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            max_neighborhood_size = max_neighborhood_size,&
            compute_influence = compute_influence,&
            save_factorization = save_factorization,&
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
    end subroutine loess_fit_robust_alloc

end module tox_loess
