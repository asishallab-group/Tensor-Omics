#include "macros.h"

!> AUTHOR_VIVIAN_BASS
module tox_loess
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, validate_in_range_int, validate_all_in_range_real, check_alloc_stat, ERR_INVALID_INPUT, ERR_ALLOC_FAIL, ERR_SIZE_MISMATCH

    ! ---- LOESS netlib externals ----
    ! ============================================================
    ! LOESS Subroutines from Netlib
    ! ============================================================
    ! These subroutines are part of the LOESS implementation provided by the Netlib library.
    ! They are used to perform various operations such as decomposition, fitting, evaluation,
    ! and robust weight computation for LOESS models. The subroutines are interfaced here
    ! to integrate the Netlib routines into the Tensor Omics project.
    !
    ! Subroutine Descriptions:
    ! - lowesd: Initializes workspace arrays and performs LOESS decomposition.
    ! - lowesb: Fits the LOESS model and computes the diagonal elements of the hat matrix.
    ! - lowese: Evaluates the LOESS model and computes the smoothed response variable array.
    ! - lowesw: Computes robust weights for LOESS using residuals.
    ! ============================================================
    interface loess_decomposition
        ! ============================================================
        ! Subroutine: lowesd
        ! ============================================================
        !> Perform LOESS decomposition.
        !| This subroutine computes the decomposition of the LOESS model.
        !| It initializes the integer and real workspace arrays based on the input parameters.
        subroutine lowesd(mode, int_workspace, int_workspace_size, real_workspace_size, real_workspace, dim, n, span, degree, nvmax, setlf)
            use, intrinsic :: iso_fortran_env, only: real64, int32
            integer(int32), intent(in) :: mode
                !! Mode of operation
            integer(int32), intent(in) :: int_workspace_size
                !! Length of the integer workspace array
            integer(int32), intent(in) :: real_workspace_size
                !! Length of the real workspace array
            integer(int32), intent(in) :: dim
                !! Dimensionality of the data
            integer(int32), intent(in) :: n
                !! Number of data points
            integer(int32), intent(in) :: degree
                !! Degree of the LOESS polynomial
            integer(int32), intent(in) :: nvmax
                !! Maximum neighborhood size
            real(real64), intent(in) :: span
                !! Smoothing parameter for LOESS
            logical, intent(in) :: setlf
                !! Save matrix factorization flag
            integer(int32), intent(inout) :: int_workspace(int_workspace_size)
                !! Integer workspace array
            real(real64), intent(inout) :: real_workspace(real_workspace_size)
                !! Real workspace array
        end subroutine lowesd
    end interface loess_decomposition

    interface loess_fitting
        ! ============================================================
        ! Subroutine: lowesb
        ! ============================================================
        !> Perform LOESS fitting.
        !| This subroutine fits the LOESS model to the data.
        !| It uses the input predictor and response variables to compute the diagonal elements of the hat matrix.
        subroutine lowesb(x, y, w, diagl, infl, int_workspace, int_workspace_size, real_workspace_size, real_workspace)
            import :: real64, int32
            integer(int32), intent(in) :: int_workspace_size
                !! Length of the integer workspace array
            integer(int32), intent(in) :: real_workspace_size
                !! Length of the real workspace array
            real(real64), intent(in) :: x(*)
                !! Predictor variable array
            real(real64), intent(in) :: y(*)
                !! Response variable array
            real(real64), intent(in) :: w(*)
                !! Weight array for data points
            real(real64), intent(out) :: diagl(*)
                !! Diagonal elements of the hat matrix
            logical, intent(in) :: infl
                !! Influence calculation flag
            integer(int32), intent(in) :: int_workspace(int_workspace_size)
                !! Integer workspace array
            real(real64), intent(inout) :: real_workspace(real_workspace_size)
                !! Real workspace array
        end subroutine lowesb
    end interface loess_fitting

    interface loess_evaluation
        ! ============================================================
        ! Subroutine: lowese
        ! ============================================================
        !> Perform LOESS evaluation.
        !| This subroutine evaluates the LOESS model.
        !| It computes the smoothed response variable array based on the input predictor variables.
        subroutine lowese(int_workspace, int_workspace_size, real_workspace_size, real_workspace, n, eval_points, fitted_values)
            import :: real64, int32
            integer(int32), intent(in) :: int_workspace_size
                !! Length of the integer workspace array
            integer(int32), intent(in) :: real_workspace_size
                !! Length of the real workspace array
            integer(int32), intent(in) :: n
                !! Number of data points
            integer(int32), intent(inout) :: int_workspace(int_workspace_size)
                !! Integer workspace array
            real(real64), intent(inout) :: real_workspace(real_workspace_size)
                !! Real workspace array
            real(real64), intent(in) :: eval_points(n, 1)
                !! x-values at which to evaluate the fitted LOESS curve (evaluation points)
            real(real64), intent(out) :: fitted_values(n)
                !! Smoothed response variable array
        end subroutine lowese
    end interface loess_evaluation

    interface loess_robust_weights
        ! ============================================================
        ! Subroutine: lowesw
        ! ============================================================
        !> Compute robust weights for LOESS.
        !| This subroutine computes the robust weights for the LOESS model.
        !| It uses the residuals to update the weights for robust fitting.
        subroutine lowesw(residuals, n, robust_weights, permutation_indices)
            use, intrinsic :: iso_fortran_env, only: real64, int32
            integer(int32), intent(in) :: n
                !! Number of data points
            real(real64), intent(in) :: residuals(n)
                !! Residuals array
            real(real64), intent(out) :: robust_weights(n)
                !! Robust weights array
            integer(int32), intent(out) :: permutation_indices(n)
                !! Permutation indices array
        end subroutine lowesw
    end interface loess_robust_weights

#define CM_MIN_REAL_WORKSPACE_SIZE 100000_int32
#define CM_MIN_INT_WORKSPACE_SIZE 10000_int32

    real(real64), parameter :: EPS_LOESS = 1.0e-12_real64

contains

    ! ============================================================
    ! Recommend workspace sizes based on Netlib exact formulas
    ! ============================================================
    !> Recommend workspace sizes based on Netlib exact formulas.
    !| Computes the required sizes for integer and real workspace arrays.
    !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
    subroutine tox_loess_required_workspace(d, nvmax, int_workspace_size, real_workspace_size, setlf)
        ! nvmax is the maximum neighborhood size (usually n_train).
        ! setlf indicates if matrix factorizations are saved.
        integer(int32), intent(in) :: d
            !! Dimensionality of the data
        integer(int32), intent(in) :: nvmax
            !! Maximum neighborhood size
        logical, intent(in) :: setlf
            !! Save matrix factorization flag
        integer(int32), intent(out) :: int_workspace_size
            !! Required size of the integer workspace array
        integer(int32), intent(out) :: real_workspace_size
            !! Required size of the real workspace array

        ! These formulas account for dimensionality (d), neighborhood size (nvmax),
        ! and whether matrix factorizations need to be saved (setlf).
        int_workspace_size = 100 + (2**d + 15) * nvmax
        if (setlf) then
            int_workspace_size = int_workspace_size + nvmax
        end if
        int_workspace_size = max(CM_MIN_INT_WORKSPACE_SIZE, int_workspace_size)

        real_workspace_size = 100 + (10*d + 20)*nvmax
        if (setlf) then 
            real_workspace_size = real_workspace_size + (d + 1) * nvmax
        end if
        real_workspace_size = max(CM_MIN_REAL_WORKSPACE_SIZE, real_workspace_size)
    end subroutine tox_loess_required_workspace

    ! ============================================================
    ! Plain LOESS fitting
    ! ============================================================
    !> Perform plain LOESS fitting.
    !| Fits a LOESS model to the data using the specified smoothing parameter.
    !| Outputs the smoothed response variable array.
    subroutine loess_fit_plain(n, x, y, w, eval_points, span, degree, nvmax, infl, setlf, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, fitted_values, ierr)
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
        integer(int32), intent(in) :: nvmax
            !! Maximum neighborhood size
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array

        real(real64), intent(in) :: x(n)
            !! Predictor variable array
        real(real64), intent(in) :: y(n)
            !! Response variable array
        real(real64), intent(in) :: w(n)
            !! Weight array for data points
        real(real64), intent(in) :: eval_points(n, 1)
            !! Evaluation points (x values at which the fitted curve is computed)
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS

        logical, intent(in) :: infl
            !! Influence calculation flag
        logical, intent(in) :: setlf
            !! Save matrix factorization flag

        integer(int32), intent(inout) :: int_workspace(int_workspace_size)
            !! Integer workspace array
        real(real64), intent(inout) :: real_workspace(real_workspace_size)
            !! Real workspace array
        real(real64), intent(inout) :: diagl(n)
            !! Diagonal elements of the hat matrix

        real(real64), intent(out) :: fitted_values(n)
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: n_eff
            !! Number of effective points

        ! Initialize error code
        call set_ok(ierr)

        ! Validate inputs
        call validate_dimension_size(n, ierr)
        call validate_in_range_real(span, ierr, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(int_workspace_size, ierr, min=CM_MIN_INT_WORKSPACE_SIZE)
        call validate_in_range_int(real_workspace_size, ierr, min=CM_MIN_REAL_WORKSPACE_SIZE)
        call validate_in_range_int(degree, ierr, min=0, max=2_int32)
        if (is_err(ierr)) return

        if (n == 1) then
            write (*, '(A)') "LOESS: Single point detected. Skipping adjustment (fitted_values = y)."
            fitted_values(1) = y(1)
            return
        end if

        ! Validate effective points
        n_eff = max(2_int32, int(ceiling(span*real(n, real64))))
        call validate_in_range_int(n_eff, ierr, min=degree + 3_int32)
        if (is_err(ierr)) return

        ! Perform the three-step LOESS fitting procedure:
        ! 1. Decomposition: Initialize workspace arrays and decompose the problem
        ! 2. Fitting: Fit the LOESS model and compute influence diagnostics
        ! 3. Evaluation: Evaluate the fitted model at eval_points to produce smoothed values
        call loess_decomposition(106, int_workspace, int_workspace_size, real_workspace_size, real_workspace, 1_int32, n, span, degree, nvmax, setlf)
        call loess_fitting(x, y, w, diagl, infl, int_workspace, int_workspace_size, real_workspace_size, real_workspace)
        call loess_evaluation(int_workspace, int_workspace_size, real_workspace_size, real_workspace, n, eval_points, fitted_values)
    end subroutine loess_fit_plain

    ! ============================================================
    ! Robust LOESS fitting
    ! ============================================================
    !> Perform robust LOESS fitting with bisquare reweighting.
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !|  - Combines original weights with robust weights (down-weights from previous iteration)
    !|  - Runs LOESS fitting with combined weights
    !|  - Computes residuals (y - fitted values)
    !|  - Updates robust weights using bisquare function (suppresses large residuals)
    !|
    subroutine loess_fit_robust(n, x, y, w, eval_points, span, degree, nvmax, infl, setlf, n_iters, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, robust_weights, combined_weights, residuals, permutation_indices, fitted_values, ierr)
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
        integer(int32), intent(in) :: nvmax
            !! Maximum neighborhood size
        integer(int32), intent(in) :: n_iters
            !! Number of robust iterations
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array

        real(real64), intent(in) :: x(n)
            !! Predictor variable array
        real(real64), intent(in) :: y(n)
            !! Response variable array
        real(real64), intent(in) :: w(n)
            !! Weight array for data points
        real(real64), intent(in) :: eval_points(n, 1)
            !! Evaluation points (x values at which the fitted curve is computed)
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS

        logical, intent(in) :: infl
            !! Influence calculation flag
        logical, intent(in) :: setlf
            !! Save matrix factorization flag

        integer(int32), intent(inout) :: int_workspace(int_workspace_size)
            !! Integer workspace array
        real(real64), intent(inout) :: real_workspace(real_workspace_size)
            !! Real workspace array
        real(real64), intent(inout) :: diagl(n)
            !! Diagonal elements of the hat matrix
        real(real64), intent(inout) :: robust_weights(n)
            !! Robust bisquare weights (updated each iteration, initialized to 1.0)
        real(real64), intent(inout) :: combined_weights(n)
            !! Combined weights: product of user weights and robust weights (w(i) * robust_weights(i))
        real(real64), intent(inout) :: residuals(n)
            !! Residuals (y - fitted_values), used to compute bisquare robust weights
        integer(int32), intent(inout) :: permutation_indices(n)
            !! Permutation indices array (from NetLib bisquare weight computation)

        real(real64), intent(out) :: fitted_values(n)
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: n_eff
            !! Number of effective points

        integer(int32) :: it, i, dim_val

        ! Initialize error code
        call set_ok(ierr)

        ! Validate inputs
        call validate_dimension_size(n, ierr)
        call validate_in_range_real(span, ierr, min=EPS_LOESS, max=1.0_real64)
        call validate_in_range_int(n_iters, ierr, min= 1_int32)
        call validate_in_range_int(int_workspace_size, ierr, min=CM_MIN_INT_WORKSPACE_SIZE)
        call validate_in_range_int(real_workspace_size, ierr, min=CM_MIN_REAL_WORKSPACE_SIZE)
        call validate_in_range_int(degree, ierr, min=0, max=2_int32)
        if (is_err(ierr)) return

        if (n == 1) then
            write (*, '(A)') "LOESS: Single point detected. Skipping adjustment (fitted_values = y)."
            fitted_values(1) = y(1)
            return
        end if

        ! Validate effective points
        n_eff = max(2_int32, int(ceiling(span*real(n, real64))))
        call validate_in_range_int(n_eff, ierr, min=degree + 3_int32)
        if (is_err(ierr)) return

        ! Initialize robust weights to 1 (no reweighting on first iteration)
        robust_weights = 1.0_real64
        dim_val = 1_int32

        ! Perform robust iterative refinement
        do it = 1, n_iters
            ! Reset workspace arrays for this iteration
            int_workspace = 0_int32
            real_workspace = 0.0_real64

            ! Combine original weights with robust weights for this iteration
            do i = 1, n
                combined_weights(i) = w(i)*robust_weights(i)
            end do

            ! Perform LOESS fitting for this robust iteration
            call loess_decomposition(106_int32, int_workspace, int_workspace_size, real_workspace_size, real_workspace, dim_val, n, span, degree, nvmax, setlf)
            call loess_fitting(x, y, combined_weights, diagl, infl, int_workspace, int_workspace_size, real_workspace_size, real_workspace)
            call loess_evaluation(int_workspace, int_workspace_size, real_workspace_size, real_workspace, n, eval_points, fitted_values)

            ! Compute residuals for robust reweighting in next iteration
            do i = 1, n
                residuals(i) = y(i) - fitted_values(i)
            end do

            ! Compute new robust weights using bisquare reweighting
            call loess_robust_weights(residuals, n, robust_weights, permutation_indices)
        end do
    end subroutine loess_fit_robust

    ! ============================================================
    ! Wrapper subroutine for LOESS fitting (plain or robust)
    ! ============================================================
    !> Wrapper subroutine for LOESS fitting (plain or robust).
    !| This subroutine selects between plain and robust LOESS fitting based on the mode.
    !| It dynamically allocates the required arrays and computes workspace sizes.
    !|
    !| Parameters:
    !| - mode: Specifies the type of LOESS fitting to perform.
    !|   - 0: Plain LOESS fitting. This mode performs a single pass of LOESS fitting without any additional weighting or iterations. It is suitable for datasets without significant outliers.
    !|   - 1: Robust LOESS fitting. This mode applies bisquare reweighting over multiple iterations to reduce the influence of outliers. The number of iterations is controlled by the `n_iters` parameter.
    subroutine loess_alloc(x, y, span, degree, fitted_values, mode, n_iters, ierr)
        use, intrinsic :: iso_fortran_env, only: real64, int32
        implicit none

        ! Input parameters
        real(real64), intent(in) :: x(:)
            !! Predictor variable array
        real(real64), intent(in) :: y(:)
            !! Response variable array
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
        integer(int32), intent(in) :: mode
            !! Mode of operation: 0 for plain, 1 for robust
        integer(int32), intent(in) :: n_iters
            !! Number of robust iterations (only used when mode = 1)

        ! Output parameters
        real(real64), intent(out) :: fitted_values(size(y))
            !! Fitted (smoothed) values of y
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local variables
        integer(int32) :: n, int_workspace_size, real_workspace_size, istat
        integer(int32), allocatable :: int_workspace(:), permutation_indices(:)
        real(real64), allocatable :: real_workspace(:), diagl(:), robust_weights(:), combined_weights(:), residuals(:), w_init(:), eval_points_mat(:, :)
        real(real64) :: range_x
        integer(int32) :: uniq_count, need_uniq
        real(real64) :: uniq_x(4)
        logical :: found
        integer(int32) :: i, j
        real(real64) :: tol

        ! Initialize variables
        n = size(y)
        call set_ok(ierr)
        call set_ok(istat)

        if (size(x) /= size(y)) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        call validate_dimension_size(n, ierr)
        if (is_err(ierr)) return

        call validate_in_range_int(mode, ierr, min=0_int32, max=1_int32)
        if (is_err(ierr)) return

        call validate_in_range_int(degree, ierr, min=0_int32, max=2_int32)
        if (is_err(ierr)) return

        call validate_in_range_real(span, ierr, min=EPS_LOESS, max=1.0_real64)
        if (is_err(ierr)) return

        if (mode == 1_int32) then
            call validate_in_range_int(n_iters, ierr, min=1_int32)
            if (is_err(ierr)) return
        end if

        if (n == 1) then
            write (*, '(A)') "LOESS: Single point detected. Skipping adjustment (fitted_values = y)."
            fitted_values(1) = y(1)
            call set_ok(ierr)
            return
        end if

        call validate_all_in_range_real(x, n, ierr)
        if (is_err(ierr)) return

        range_x = maxval(x) - minval(x)
        if (range_x <= EPS_LOESS) then
            write (*, '(A, E12.4)') "LOESS: Range of x is too small (<= EPS). Skipping adjustment. Range =", range_x
            fitted_values = y
            call set_ok(ierr)
            return
        end if

        need_uniq = min(4_int32, degree + 2_int32)
        uniq_count = 0_int32
        uniq_x = 0.0_real64

        do i = 1, n
            ! tolerance to compare floats
            tol = EPS_LOESS*max(1.0_real64, abs(x(i)))

            found = .false.
            do j = 1, uniq_count
                if (abs(x(i) - uniq_x(j)) <= tol) then
                    found = .true.
                    exit
                end if
            end do

            if (.not. found) then
                uniq_count = uniq_count + 1_int32
                uniq_x(uniq_count) = x(i)
                if (uniq_count >= need_uniq) exit
            end if
        end do

        if (uniq_count < need_uniq) then
            ! a lot of same values to adjust
            write (*, '(A, I2, A, I2, A)') "LOESS: Insufficient unique points (Found ", uniq_count, &
                " but need ", need_uniq, "). Using identity mapping."
            fitted_values = y
            return
        end if

        call tox_loess_required_workspace(1_int32, n, int_workspace_size, real_workspace_size, .false.)

        ! Allocate workspace arrays
        allocate (int_workspace(int_workspace_size), real_workspace(real_workspace_size), diagl(n), w_init(n), eval_points_mat(n, 1), stat=istat)
        call check_alloc_stat(istat, ierr)
        if (is_err(ierr)) return

        ! Initialize arrays
        int_workspace = 0_int32
        w_init = 1.0_real64  ! Uniform weights initially
        diagl = 0.0_real64
        eval_points_mat(:, 1) = x  ! Evaluate at training x-values to obtain in-sample fitted values
        real_workspace = 0.0_real64

        ! Allocate additional arrays for robust mode
        if (mode == 1_int32) then
            allocate (robust_weights(n), combined_weights(n), residuals(n), permutation_indices(n), stat=istat)
            call check_alloc_stat(istat, ierr)
            if (is_err(ierr)) then
                deallocate (int_workspace, real_workspace, diagl, w_init, eval_points_mat)
                return
            end if
        end if

        ! Call the appropriate LOESS fitting subroutine based on mode
        ! Mode 0: Plain fitting (single pass, no robust iterations)
        ! Mode 1: Robust fitting (multiple iterations with reweighting)
        if (mode == 0_int32) then
            call loess_fit_plain(n, x, y, w_init, eval_points_mat, span, degree, n, &
                                 .false., .false., int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, fitted_values, ierr)
        else
            call loess_fit_robust(n, x, y, w_init, eval_points_mat, span, degree, n, &
                                  .false., .false., n_iters, int_workspace, int_workspace_size, real_workspace, real_workspace_size, &
                                  diagl, robust_weights, combined_weights, residuals, permutation_indices, fitted_values, ierr)
        end if

        ! Deallocate workspace arrays
        deallocate (int_workspace, real_workspace, diagl, w_init, eval_points_mat)
        if (allocated(robust_weights)) deallocate (robust_weights, combined_weights, residuals, permutation_indices)

    end subroutine loess_alloc

end module tox_loess

!> C interface for the LOESS wrapper subroutine.
!| This subroutine selects between plain and robust LOESS fitting based on the mode.
!| It dynamically allocates the required arrays and computes workspace sizes.
!|
!| Parameters:
!| - mode: Specifies the type of LOESS fitting to perform.
!|   - 0: Plain LOESS fitting. This mode performs a single pass of LOESS fitting without any additional weighting or iterations. It is suitable for datasets without significant outliers.
!|   - 1: Robust LOESS fitting. This mode applies bisquare reweighting over multiple iterations to reduce the influence of outliers. The number of iterations is controlled by the `n_iters` parameter.
subroutine tox_loess_c(x, y, n, span, degree, fitted_values, mode, n_iters, ierr) bind(C, name="tox_loess_c")
    use tox_loess, only: loess_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    ! Input parameters
    integer(c_int), intent(in), target :: n
        !! Number of data points
    real(c_double), dimension(n), intent(in), target :: x
        !! Predictor variable array
    real(c_double), dimension(n), intent(in), target :: y
        !! Response variable array
    real(c_double), intent(in), target :: span
        !! Smoothing parameter for LOESS
    integer(c_int), intent(in), target :: degree
        !! Degree of the LOESS polynomial
    integer(c_int), intent(in), target :: mode
        !! Mode of operation: 0 for plain, 1 for robust
    integer(c_int), intent(in), target :: n_iters
        !! Number of robust iterations (only used when mode = 1)

    ! Output parameters
    real(c_double), dimension(n), intent(out), target :: fitted_values
        !! Fitted (smoothed) values of y
    integer(c_int), intent(out), target :: ierr
        !! Error code

    ! Null validation macros
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n)
    M_CHECK_NON_NULL(x)
    M_CHECK_NON_NULL(y)
    M_CHECK_NON_NULL(span)
    M_CHECK_NON_NULL(degree)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(n_iters)
    M_CHECK_NON_NULL(fitted_values)

    ! Call the Fortran subroutine
    call loess_alloc(x, y, span, degree, fitted_values, mode, n_iters, ierr)
end subroutine tox_loess_c

!> Perform plain LOESS fitting.
!| Fits a LOESS model to the data using the specified smoothing parameter.
!| Outputs the smoothed response variable array.
subroutine loess_fit_plain_c(n, x, y, w, eval_points, span, degree, nvmax, infl, setlf, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, fitted_values, ierr) bind(C, name="loess_fit_plain_c")
    use tox_loess, only: loess_fit_plain
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    ! Input parameters
    integer(c_int), intent(in), target :: n
        !! Total number of data points
    real(c_double), dimension(n), intent(in), target :: x
        !! Predictor variable array
    real(c_double), dimension(n), intent(in), target :: y
        !! Response variable array
    real(c_double), dimension(n), intent(in), target :: w
        !! Weight array for data points
    real(c_double), dimension(n, 1), intent(in), target :: eval_points
        !! Evaluation points (x values at which the fitted curve is computed)
    real(c_double), intent(in), target :: span
        !! Smoothing parameter for LOESS
    integer(c_int), intent(in), target :: degree
        !! Degree of the LOESS polynomial
    integer(c_int), intent(in), target :: nvmax
        !! Maximum neighborhood size
    integer(c_int), intent(in), target :: infl
        !! Influence calculation flag (0 for false, non-zero for true)
    integer(c_int), intent(in), target :: setlf
        !! Save matrix factorization flag (0 for false, non-zero for true)
    integer(c_int), dimension(*), intent(inout), target :: int_workspace
        !! Integer workspace array
    integer(c_int), intent(in), target :: int_workspace_size
        !! Required size of the integer workspace array
    real(c_double), dimension(*), intent(inout), target :: real_workspace
        !! Real workspace array
    integer(c_int), intent(in), target :: real_workspace_size
        !! Required size of the real workspace array

    ! Output parameters
    real(c_double), dimension(n), intent(out), target :: diagl
        !! Diagonal elements of the hat matrix
    real(c_double), dimension(n), intent(out), target :: fitted_values
        !! Fitted (smoothed) values of y at the evaluation points
    integer(c_int), intent(out), target :: ierr
        !! Error code

    ! Null validation macros
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

    ! Call the Fortran subroutine
    call loess_fit_plain(n, x, y, w, eval_points, span, degree, nvmax, infl /= 0, setlf /= 0, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, fitted_values, ierr)
end subroutine loess_fit_plain_c

!> Perform robust LOESS fitting with bisquare reweighting.
!| Fits a LOESS model to the data using robust iterations to handle outliers.
!| Outputs the smoothed response variable array.
subroutine loess_fit_robust_c(n, x, y, w, eval_points, span, degree, nvmax, infl, setlf, n_iters, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, robust_weights, combined_weights, residuals, permutation_indices, fitted_values, ierr) bind(C, name="loess_fit_robust_c")
    use tox_loess, only: loess_fit_robust
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    ! Input parameters
    integer(c_int), intent(in), target :: n
        !! Number of data points
    real(c_double), dimension(n), intent(in), target :: x
        !! Predictor variable array
    real(c_double), dimension(n), intent(in), target :: y
        !! Response variable array
    real(c_double), dimension(n), intent(in), target :: w
        !! Weight array for data points
    real(c_double), dimension(n, 1), intent(in), target :: eval_points
        !! Evaluation points (x values at which the fitted curve is computed)
    real(c_double), intent(in), target :: span
        !! Smoothing parameter for LOESS
    integer(c_int), intent(in), target :: degree
        !! Degree of the LOESS polynomial
    integer(c_int), intent(in), target :: nvmax
        !! Maximum neighborhood size
    integer(c_int), intent(in), target :: infl
        !! Influence calculation flag (0 for false, non-zero for true)
    integer(c_int), intent(in), target :: setlf
        !! Save matrix factorization flag (0 for false, non-zero for true)
    integer(c_int), intent(in), target :: n_iters
        !! Number of robust iterations
    integer(c_int), dimension(*), intent(inout), target :: int_workspace
        !! Integer workspace array
    integer(c_int), intent(in), target :: int_workspace_size
        !! Required size of the integer workspace array
    real(c_double), dimension(*), intent(inout), target :: real_workspace
        !! Real workspace array
    integer(c_int), intent(in), target :: real_workspace_size
        !! Required size of the real workspace array

    ! Output parameters
    real(c_double), dimension(n), intent(out), target :: diagl
        !! Diagonal elements of the hat matrix
    real(c_double), dimension(n), intent(out), target :: robust_weights
        !! Robust bisquare weights (updated each iteration)
    real(c_double), dimension(n), intent(out), target :: combined_weights
        !! Combined weights: product of user weights and robust weights
    real(c_double), dimension(n), intent(out), target :: residuals
        !! Residuals (y - fitted_values), used to compute robust weights
    integer(c_int), dimension(n), intent(out), target :: permutation_indices
        !! Permutation indices array (from NetLib bisquare weight computation)
    real(c_double), dimension(n), intent(out), target :: fitted_values
        !! Fitted (smoothed) values of y at the evaluation points
    integer(c_int), intent(out), target :: ierr
        !! Error code

    ! Null validation macros
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

    ! Call the Fortran subroutine
    call loess_fit_robust(n, x, y, w, eval_points, span, degree, nvmax, infl /= 0, setlf /= 0, n_iters, int_workspace, int_workspace_size, real_workspace, real_workspace_size, diagl, robust_weights, combined_weights, residuals, permutation_indices, fitted_values, ierr)
end subroutine loess_fit_robust_c

!! Wrapper for recommending workspace sizes based on Netlib exact formulas.
!! This subroutine is designed to be called from C code and computes the required sizes
!! for integer and real workspace arrays. These sizes depend on the dimensionality of the data,
!! the maximum neighborhood size, and whether matrix factorizations are saved.
subroutine tox_loess_required_workspace_c(d, nvmax, int_workspace_size, real_workspace_size, setlf) bind(C, name="tox_loess_required_workspace_c")
    use tox_loess, only: tox_loess_required_workspace
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none

    ! Input parameters
    integer(c_int), intent(in), target :: d
        !! Dimensionality of the data
    integer(c_int), intent(in), target :: nvmax
        !! Maximum neighborhood size
    integer(c_int), intent(in), target :: setlf
        !! Save matrix factorization flag (0 for false, non-zero for true)

    ! Output parameters
    integer(c_int), intent(out) :: int_workspace_size
        !! Required size of the integer workspace array
    integer(c_int), intent(out) :: real_workspace_size
        !! Required size of the real workspace array

    ! Call the Fortran subroutine
    call tox_loess_required_workspace(d, nvmax, int_workspace_size, real_workspace_size, setlf /= 0)
end subroutine tox_loess_required_workspace_c
