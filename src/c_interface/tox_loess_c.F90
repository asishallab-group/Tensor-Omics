#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_loess(module)]]
module tox_loess_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[tox_loess(module):tox_loess_required_workspace(subroutine)]]
    !| Recommend workspace sizes based on Netlib exact formulas.
    !| Computes the required sizes for integer and real workspace arrays.
    !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
    subroutine tox_loess_required_workspace_c(&
            d,&
            nvmax,&
            int_workspace_size,&
            real_workspace_size,&
            setlf,&
            ierr&
            ) bind(C, name="tox_loess_required_workspace_c")
        use tox_loess, only: tox_loess_required_workspace
        use tox_loess

        integer(c_int), intent(in), target :: d
            !! Dimensionality of the data
        integer(c_int), intent(in), target :: nvmax
            !! Maximum neighborhood size
        integer(c_int), intent(out), target :: int_workspace_size
            !! Required size of the integer workspace array
        integer(c_int), intent(out), target :: real_workspace_size
            !! Required size of the real workspace array
        integer(c_int), intent(in), target :: setlf
            !! Save matrix factorization flag
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: setlf_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(d)
        M_CHECK_NON_NULL(nvmax)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(setlf)

        call c_int_as_logical(setlf, setlf_f)
        call tox_loess_required_workspace(&
            d = d,&
            nvmax = nvmax,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            setlf = setlf_f&
        )

    end subroutine tox_loess_required_workspace_c

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_plain(subroutine)]]
    !| Perform plain LOESS fitting.
    !| Fits a LOESS model to the data using the specified smoothing parameter.
    !| Outputs the smoothed response variable array.
    subroutine loess_fit_plain_c(&
            n,&
            x,&
            y,&
            w,&
            eval_points,&
            span,&
            degree,&
            nvmax,&
            infl,&
            setlf,&
            int_workspace,&
            int_workspace_size,&
            real_workspace,&
            real_workspace_size,&
            diagl,&
            fitted_values,&
            ierr&
            ) bind(C, name="loess_fit_plain_c")
        use tox_loess, only: loess_fit_plain
        use tox_loess
        integer(c_int), intent(in), target :: n
            !! Total number of data points
        integer(c_int), intent(in), target :: int_workspace_size
            !! Required size of the integer workspace array
        integer(c_int), intent(in), target :: real_workspace_size
            !! Required size of the real workspace array
        real(c_double), intent(in), dimension(n), target :: x
            !! Predictor variable array
        real(c_double), intent(in), dimension(n), target :: y
            !! Response variable array
        real(c_double), intent(in), dimension(n), target :: w
            !! Weight array for data points
        real(c_double), intent(in), dimension(n, 1), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
        integer(c_int), intent(in), target :: nvmax
            !! Maximum neighborhood size
        integer(c_int), intent(in), target :: infl
            !! Influence calculation flag
        integer(c_int), intent(in), target :: setlf
            !! Save matrix factorization flag
        integer(c_int), intent(inout), dimension(int_workspace_size), target :: int_workspace
            !! Integer workspace array
        real(c_double), intent(inout), dimension(real_workspace_size), target :: real_workspace
            !! Real workspace array
        real(c_double), intent(inout), dimension(n), target :: diagl
            !! Diagonal elements of the hat matrix
        real(c_double), intent(out), dimension(n), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: infl_f
        logical :: setlf_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(x)
        M_CHECK_NON_NULL(y)
        M_CHECK_NON_NULL(w)
        M_CHECK_NON_NULL(eval_points)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(nvmax)
        M_CHECK_NON_NULL(infl)
        M_CHECK_NON_NULL(setlf)
        M_CHECK_NON_NULL(int_workspace)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(diagl)
        M_CHECK_NON_NULL(fitted_values)

        call c_int_as_logical(infl, infl_f)
        call c_int_as_logical(setlf, setlf_f)
        call loess_fit_plain(&
            n = n,&
            x = x,&
            y = y,&
            w = w,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            nvmax = nvmax,&
            infl = infl_f,&
            setlf = setlf_f,&
            int_workspace = int_workspace,&
            int_workspace_size = int_workspace_size,&
            real_workspace = real_workspace,&
            real_workspace_size = real_workspace_size,&
            diagl = diagl,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )

    end subroutine loess_fit_plain_c

    !> summary: C-wrapper for [[tox_loess(module):loess_fit_robust(subroutine)]]
    !| Perform robust LOESS fitting with bisquare reweighting.
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
            w,&
            eval_points,&
            span,&
            degree,&
            nvmax,&
            infl,&
            setlf,&
            n_iters,&
            int_workspace,&
            int_workspace_size,&
            real_workspace,&
            real_workspace_size,&
            diagl,&
            robust_weights,&
            combined_weights,&
            residuals,&
            permutation_indices,&
            fitted_values,&
            ierr&
            ) bind(C, name="loess_fit_robust_c")
        use tox_loess, only: loess_fit_robust
        use tox_loess
        integer(c_int), intent(in), target :: n
            !! Total number of data points
        integer(c_int), intent(in), target :: int_workspace_size
            !! Required size of the integer workspace array
        integer(c_int), intent(in), target :: real_workspace_size
            !! Required size of the real workspace array
        real(c_double), intent(in), dimension(n), target :: x
            !! Predictor variable array
        real(c_double), intent(in), dimension(n), target :: y
            !! Response variable array
        real(c_double), intent(in), dimension(n), target :: w
            !! Weight array for data points
        real(c_double), intent(in), dimension(n, 1), target :: eval_points
            !! Evaluation points (x values at which the fitted curve is computed)
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
        integer(c_int), intent(in), target :: nvmax
            !! Maximum neighborhood size
        integer(c_int), intent(in), target :: infl
            !! Influence calculation flag
        integer(c_int), intent(in), target :: setlf
            !! Save matrix factorization flag
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations
        integer(c_int), intent(inout), dimension(int_workspace_size), target :: int_workspace
            !! Integer workspace array
        real(c_double), intent(inout), dimension(real_workspace_size), target :: real_workspace
            !! Real workspace array
        real(c_double), intent(inout), dimension(n), target :: diagl
            !! Diagonal elements of the hat matrix
        real(c_double), intent(inout), dimension(n), target :: robust_weights
            !! Robust bisquare weights (updated each iteration, initialized to 1.0)
        real(c_double), intent(inout), dimension(n), target :: combined_weights
            !! Combined weights: product of user weights and robust weights (w(i) * robust_weights(i))
        real(c_double), intent(inout), dimension(n), target :: residuals
            !! Residuals (y - fitted_values), used to compute bisquare robust weights
        integer(c_int), intent(inout), dimension(n), target :: permutation_indices
            !! Permutation indices array (from NetLib bisquare weight computation)
        real(c_double), intent(out), dimension(n), target :: fitted_values
            !! Fitted (smoothed) values of y at the evaluation points
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical :: infl_f
        logical :: setlf_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(x)
        M_CHECK_NON_NULL(y)
        M_CHECK_NON_NULL(w)
        M_CHECK_NON_NULL(eval_points)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(nvmax)
        M_CHECK_NON_NULL(infl)
        M_CHECK_NON_NULL(setlf)
        M_CHECK_NON_NULL(n_iters)
        M_CHECK_NON_NULL(int_workspace)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(diagl)
        M_CHECK_NON_NULL(robust_weights)
        M_CHECK_NON_NULL(combined_weights)
        M_CHECK_NON_NULL(residuals)
        M_CHECK_NON_NULL(permutation_indices)
        M_CHECK_NON_NULL(fitted_values)

        call c_int_as_logical(infl, infl_f)
        call c_int_as_logical(setlf, setlf_f)
        call loess_fit_robust(&
            n = n,&
            x = x,&
            y = y,&
            w = w,&
            eval_points = eval_points,&
            span = span,&
            degree = degree,&
            nvmax = nvmax,&
            infl = infl_f,&
            setlf = setlf_f,&
            n_iters = n_iters,&
            int_workspace = int_workspace,&
            int_workspace_size = int_workspace_size,&
            real_workspace = real_workspace,&
            real_workspace_size = real_workspace_size,&
            diagl = diagl,&
            robust_weights = robust_weights,&
            combined_weights = combined_weights,&
            residuals = residuals,&
            permutation_indices = permutation_indices,&
            fitted_values = fitted_values,&
            ierr = ierr&
        )

    end subroutine loess_fit_robust_c

    !> summary: C-wrapper for [[tox_loess(module):loess_alloc(subroutine)]]
    !| Wrapper subroutine for LOESS fitting (plain or robust).
    !| This subroutine selects between plain and robust LOESS fitting based on the mode.
    !| It dynamically allocates the required arrays and computes workspace sizes.
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
        use tox_loess
        integer(c_int), intent(in), target :: n_x_elements
            !! Size of the 1. dimension/extent of `x`
        integer(c_int), intent(in), target :: n_y_elements
            !! Size of the 1. dimension/extent of `y`
        real(c_double), intent(in), dimension(n_x_elements), target :: x
            !! Predictor variable array
        real(c_double), intent(in), dimension(n_y_elements), target :: y
            !! Response variable array
        real(c_double), intent(in), target :: span
            !! Smoothing parameter for LOESS
        integer(c_int), intent(in), target :: degree
            !! Degree of the LOESS polynomial
        real(c_double), intent(out), dimension(size(y)), target :: fitted_values
            !! Fitted (smoothed) values of y
        character(len=1, kind=c_char), intent(in), dimension(6), target :: mode
            !! Mode of operation
            !! 
            !! |      Mode      |  Value   |
            !! |----------------|----------|
            !! | robust fitting | "robust" |
            !! | plain fitting  | "plain"  |
        integer(c_int), intent(in), target :: n_iters
            !! Number of robust iterations (only used when mode = 1)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: mode_f
        integer(int32) :: mode_int_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(x)
        M_CHECK_NON_NULL(n_x_elements)
        M_CHECK_NON_NULL(y)
        M_CHECK_NON_NULL(n_y_elements)
        M_CHECK_NON_NULL(span)
        M_CHECK_NON_NULL(degree)
        M_CHECK_NON_NULL(fitted_values)
        M_CHECK_NON_NULL(mode)
        M_CHECK_NON_NULL(n_iters)
        call c_char_1d_as_string(mode, mode_f, ierr)
        if (is_err(ierr)) return

        select case (mode_f)
            case ("robust")
                    mode_int_f = MODE_ROBUST
            case ("plain")
                    mode_int_f = MODE_PLAIN
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
                return
        end select

        call loess_alloc(&
            x = x,&
            y = y,&
            span = span,&
            degree = degree,&
            fitted_values = fitted_values,&
            mode = mode_int_f,&
            n_iters = n_iters,&
            ierr = ierr&
        )

    end subroutine loess_c

end module tox_loess_c
#endif