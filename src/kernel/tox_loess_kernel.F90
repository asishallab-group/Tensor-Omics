#include <src/macros.h>

!> Kernels for LOESS (netlib `dloess`/`lowesd` family) local polynomial regression smoothing.
!| The generator turns `loess_fit_plain_kernel` / `loess_fit_robust_kernel` into the expert fitting
!| wrappers `loess_fit_plain` / `loess_fit_robust`, each with an allocating sibling, in module
!| `tox_loess`. Both name `loess_degenerate_fit` as their prologue, so data too degenerate to fit is
!| answered there and the netlib call is skipped -- the kernels themselves fit, and assume they were
!| given something fittable. There is no combined entry point that dispatches on a mode: a caller
!| chooses the plain or the robust routine, and supplies the weights and evaluation points it wants.
!| The netlib
!| interface blocks, the mode/iteration constants, `EPS_LOESS` and the workspace-sizing routine
!| `tox_loess_required_workspace` live here for callers (and the generated wrappers) to use; the
!| netlib routines themselves are not re-documented beyond their calling convention.
module tox_loess_kernel
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, validate_in_range_int, validate_all_in_range_real, check_io_stat, ERR_INVALID_INPUT, ERR_ALLOC_FAIL, ERR_SIZE_MISMATCH
    use f42_utils, only: is_close
    M_IMPLICIT_NONE

#define CM_MODE_PLAIN 0_int32
#define CM_MODE_ROBUST 1_int32
#define CM_DEFAULT_LOESS_ITERS 3_int32

    integer(int32), parameter, public :: DEFAULT_LOESS_ITERS = CM_DEFAULT_LOESS_ITERS
        !! Robust iterations used when the caller does not say, matching the count
        !! [[tox_get_outliers(module)]] and [[tox_normalization(module)]] use
    integer(int32), parameter, public :: MODE_PLAIN = CM_MODE_PLAIN
        !! Mode code selecting a plain LOESS fit, for callers that carry the choice as a value
    integer(int32), parameter, public :: MODE_ROBUST = CM_MODE_ROBUST
        !! Mode code selecting a robust LOESS fit, for callers that carry the choice as a value

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
        subroutine lowesd(iv19_code, int_workspace, int_workspace_size, real_workspace_size, real_workspace, n_dim, n, span, degree, max_neighborhood_size, save_factorization)
            use, intrinsic :: iso_fortran_env, only: real64, int32
            integer(int32), intent(in) :: iv19_code
                !! Packed netlib `iv(19)` model-selection code (family/surface/statistics digits); see
                !! [[tox_loess(module):loess_fit_plain(subroutine)]] for the fixed value used by this codebase.
            integer(int32), intent(in) :: int_workspace_size
                !! Length of the integer workspace array
            integer(int32), intent(in) :: real_workspace_size
                !! Length of the real workspace array
            integer(int32), intent(in) :: n_dim
                !! Dimensionality of the data
            integer(int32), intent(in) :: n
                !! Number of data points
            integer(int32), intent(in) :: degree
                !! Degree of the LOESS polynomial
            integer(int32), intent(in) :: max_neighborhood_size
                !! Maximum neighborhood size
            real(real64), intent(in) :: span
                !! Smoothing parameter for LOESS
            logical, intent(in) :: save_factorization
                !! Save matrix factorization flag
            integer(int32), intent(out) :: int_workspace(int_workspace_size)
                !! Integer workspace array, laid out and populated here; nothing in it is read
            real(real64), intent(out) :: real_workspace(real_workspace_size)
                !! Real workspace array, whose first entries are set here; nothing in it is read
        end subroutine lowesd
    end interface loess_decomposition

    interface loess_fitting
        ! ============================================================
        ! Subroutine: lowesb
        ! ============================================================
        !> Perform LOESS fitting.
        !| This subroutine fits the LOESS model to the data.
        !| It uses the input predictor and response variables to compute the diagonal elements of the hat matrix.
        subroutine lowesb(x, y, weights, hat_diag, compute_influence, int_workspace, int_workspace_size, real_workspace_size, real_workspace)
            import :: real64, int32
            integer(int32), intent(in) :: int_workspace_size
                !! Length of the integer workspace array
            integer(int32), intent(in) :: real_workspace_size
                !! Length of the real workspace array
            real(real64), intent(in) :: x(*)
                !! Predictor variable array
            real(real64), intent(in) :: y(*)
                !! Response variable array
            real(real64), intent(in) :: weights(*)
                !! Weight array for data points
            real(real64), intent(out) :: hat_diag(*)
                !! Diagonal elements of the hat matrix
            logical, intent(in) :: compute_influence
                !! Influence calculation flag
            integer(int32), intent(inout) :: int_workspace(int_workspace_size)
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
            integer(int32), intent(in) :: int_workspace(int_workspace_size)
                !! Integer workspace array, read as laid out by the decomposition and the fit
            real(real64), intent(in) :: real_workspace(real_workspace_size)
                !! Real workspace array, read as left by the fit
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
        !! Minimum allowed LOESS `span` and a general-purpose small epsilon used across this module
        !! (and by callers such as [[tox_get_outliers(module):compute_family_scaling(subroutine)]]) to
        !! guard against division/log by (near-)zero.

contains

    !> summary: Decides whether the data is too degenerate to fit, and answers it directly if so
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The prologue of both LOESS fitting kernels. A single point, an `x` range below
    !| `EPS_LOESS`, or fewer distinct `x` values than the polynomial degree needs cannot
    !| produce a meaningful fit; rather than hand netlib an input it cannot answer, the fitted
    !| values are the observations themselves and the call reports `handled`, so the kernel is
    !| skipped. This is policy, not validation, which is why it lives here and not in the
    !| kernels: they fit, and assume they were given something fittable.
    pure subroutine loess_degenerate_fit(n, x, y, degree, fitted_values, handled, ierr)
        integer(int32), intent(in) :: n
            !! Total number of data points
        real(real64), dimension(n), intent(in) :: x
            !! Predictor variable array
        real(real64), dimension(n), intent(in) :: y
            !! Response variable array
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
        real(real64), dimension(n), intent(out) :: fitted_values
            !! Fitted values, written only when the fit was answered here
        logical, intent(out) :: handled
            !! `.true.` when the data was degenerate and `fitted_values` already holds the answer
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: tol
        real(real64) :: uniq_x(4)
        integer(int32) :: uniq_count, need_uniq, i, j
        logical :: found

        call set_ok(ierr)
        handled = .true.

        ! A single point carries no trend to smooth
        if (n == 1) then
            fitted_values(1) = y(1)
            return
        end if

        ! Near-constant x: every neighbourhood is the whole sample, so the fit is the data.
        ! Closeness is judged on LOESS's own terms -- EPS_LOESS, the smoothing floor -- rather
        ! than on what the arithmetic can resolve.
        if (is_close(maxval(x), minval(x), EPS_LOESS)) then
            fitted_values = y
            return
        end if

        ! A degree-d local polynomial needs at least d+1 distinct x-support points to be well-defined;
        ! capped at 4 (the fixed size of uniq_x below) since degree never exceeds 2 in this codebase.
        need_uniq = min(4_int32, degree + 2_int32)
        uniq_count = 0_int32
        uniq_x = 0.0_real64

        do i = 1, n
            ! Scale-relative tolerance (not a fixed epsilon) so comparisons remain meaningful
            ! for x-values far from zero.
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
            fitted_values = y
            return
        end if

        handled = .false.
    end subroutine loess_degenerate_fit

    ! ============================================================
    ! Recommend workspace sizes based on Netlib exact formulas
    ! ============================================================
    !> M_EXPORT_C
    !| summary: Recommend workspace sizes based on Netlib exact formulas
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Computes the required sizes for integer and real workspace arrays.
    !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
    subroutine tox_loess_required_workspace(n_dim, max_neighborhood_size, int_workspace_size, real_workspace_size, save_factorization)
        ! max_neighborhood_size is the maximum neighborhood size (usually n_train).
        ! save_factorization indicates if matrix factorizations are saved.
        integer(int32), intent(in) :: n_dim
            !! Dimensionality of the data
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        logical, intent(in) :: save_factorization
            !! Save matrix factorization flag
        integer(int32), intent(out) :: int_workspace_size
            !! Required size of the integer workspace array
        integer(int32), intent(out) :: real_workspace_size
            !! Required size of the real workspace array

        ! These formulas account for dimensionality (n_dim), neighborhood size (max_neighborhood_size),
        ! and whether matrix factorizations need to be saved (save_factorization).
        int_workspace_size = 100 + (2**n_dim + 15) * max_neighborhood_size
        if (save_factorization) then
            int_workspace_size = int_workspace_size + max_neighborhood_size
        end if
        int_workspace_size = max(CM_MIN_INT_WORKSPACE_SIZE, int_workspace_size)

        real_workspace_size = 100 + (10*n_dim + 20)*max_neighborhood_size
        if (save_factorization) then
            real_workspace_size = real_workspace_size + (n_dim + 1) * max_neighborhood_size
        end if
        real_workspace_size = max(CM_MIN_REAL_WORKSPACE_SIZE, real_workspace_size)
    end subroutine tox_loess_required_workspace

    ! ============================================================
    ! Plain LOESS fitting
    ! ============================================================
    !> summary: Perform plain LOESS fitting
    !| AUTHOR_FRANZ_ERIC_SILL
    !| DM_PROLOGUE(loess_degenerate_fit, tox_loess_kernel, BOTH)
    !| Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    !| response array.
    subroutine loess_fit_plain_kernel(n, x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_hat_diag, fitted_values, ierr)
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! DM_MIN(0_int32)
            !! DM_MAX(2_int32)
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_kernel, AUTO)
            !! DM_MIN(CM_MIN_INT_WORKSPACE_SIZE)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! |       n_dim       |   1_int32   |
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_kernel, AUTO)
            !! DM_MIN(CM_MIN_REAL_WORKSPACE_SIZE)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! |       n_dim       |   1_int32   |
        real(real64), intent(in) :: x(n)
            !! Predictor variable array
        real(real64), intent(in) :: y(n)
            !! Response variable array
        real(real64), intent(in) :: weights(n)
            !! Weight array for data points
        real(real64), intent(in) :: eval_points(n, 1)
            !! Evaluation points (x values at which the fitted curve is computed)
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! DM_MIN(EPS_LOESS)
            !! DM_MAX(1.0_real64)

        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! DM_DEFAULT(.false.)
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! DM_DEFAULT(.false.)

        integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
            !! Integer workspace array
        real(real64), intent(out) :: tmp_real_workspace(real_workspace_size)
            !! Real workspace array
        real(real64), intent(out) :: tmp_hat_diag(n)
            !! Diagonal elements of the hat matrix

        real(real64), intent(out) :: fitted_values(n)
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: neighborhood_size
            !! Size of the local neighborhood used at the current span (points per local fit)

        logical :: actual_compute_influence, actual_save_factorization

        call set_ok(ierr)

        M_DEFAULT_VAL(compute_influence, actual_compute_influence, .false.)
        M_DEFAULT_VAL(save_factorization, actual_save_factorization, .false.)

        ! `span` is a fraction of `n` points included in each local neighborhood; neighborhood_size is that
        ! neighborhood size. Require at least degree+3 points per neighborhood so the local
        ! polynomial fit (degree+1 coefficients) is over-determined rather than exactly/under-determined.
        neighborhood_size = max(2_int32, int(ceiling(span*real(n, real64))))
        call validate_in_range_int(neighborhood_size, ierr, min=degree + 3_int32)
        if (is_err(ierr)) return

        ! Perform the three-step LOESS fitting procedure:
        ! 1. Decomposition: Initialize workspace arrays and decompose the problem
        ! 2. Fitting: Fit the LOESS model and compute influence diagnostics
        ! 3. Evaluation: Evaluate the fitted model at eval_points to produce smoothed values
        ! `106` is netlib's packed `iv(19)` model-selection code (family/surface/statistics digits);
        ! it selects the netlib default (Gaussian family, direct surface, exact statistics) and is not
        ! meant to be tuned per call.
        call loess_decomposition(106, tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace, 1_int32, n, span, degree, max_neighborhood_size, actual_save_factorization)
        call loess_fitting(x, y, weights, tmp_hat_diag, actual_compute_influence, tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace)
        call loess_evaluation(tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace, n, eval_points, fitted_values)
    end subroutine loess_fit_plain_kernel

    ! ============================================================
    ! Robust LOESS fitting
    ! ============================================================
    !> summary: Perform robust LOESS fitting with bisquare reweighting
    !| AUTHOR_FRANZ_ERIC_SILL
    !| DM_PROLOGUE(loess_degenerate_fit, tox_loess_kernel, BOTH)
    !| Fits a LOESS model to the data using robust iterations to handle outliers.
    !| The robust fitting process iterates n_iters times, each iteration:
    !|  - Combines original weights with robust weights (down-weights from previous iteration)
    !|  - Runs LOESS fitting with combined weights
    !|  - Computes residuals (y - fitted values)
    !|  - Updates robust weights using bisquare function (suppresses large residuals)
    !|
    subroutine loess_fit_robust_kernel(n, x, y, weights, eval_points, span, degree, max_neighborhood_size, compute_influence, save_factorization, n_iters, tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_hat_diag, tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, fitted_values, ierr)
        integer(int32), intent(in) :: n
            !! Total number of data points
        integer(int32), intent(in) :: degree
            !! Degree of the LOESS polynomial
            !! DM_MIN(0_int32)
            !! DM_MAX(2_int32)
        integer(int32), intent(in) :: max_neighborhood_size
            !! Maximum neighborhood size
        integer(int32), intent(in), optional :: n_iters
            !! Number of robust iterations
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_DEFAULT_LOESS_ITERS)
        integer(int32), intent(in) :: int_workspace_size
            !! Required size of the integer workspace array
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_kernel, AUTO)
            !! DM_MIN(CM_MIN_INT_WORKSPACE_SIZE)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! |       n_dim       |   1_int32   |
        integer(int32), intent(in) :: real_workspace_size
            !! Required size of the real workspace array
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_kernel, AUTO)
            !! DM_MIN(CM_MIN_REAL_WORKSPACE_SIZE)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! |       n_dim       |   1_int32   |

        real(real64), intent(in) :: x(n)
            !! Predictor variable array
        real(real64), intent(in) :: y(n)
            !! Response variable array
        real(real64), intent(in) :: weights(n)
            !! Weight array for data points
        real(real64), intent(in) :: eval_points(n, 1)
            !! Evaluation points (x values at which the fitted curve is computed)
        real(real64), intent(in) :: span
            !! Smoothing parameter for LOESS
            !! DM_MIN(EPS_LOESS)
            !! DM_MAX(1.0_real64)

        logical, intent(in), optional :: compute_influence
            !! Influence calculation flag
            !! DM_DEFAULT(.false.)
        logical, intent(in), optional :: save_factorization
            !! Save matrix factorization flag
            !! DM_DEFAULT(.false.)

        integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
            !! Integer workspace array
        real(real64), intent(out) :: tmp_real_workspace(real_workspace_size)
            !! Real workspace array
        real(real64), intent(out) :: tmp_hat_diag(n)
            !! Diagonal elements of the hat matrix
        real(real64), intent(out) :: tmp_robust_weights(n)
            !! Robust bisquare weights (updated each iteration, initialized to 1.0)
        real(real64), intent(out) :: tmp_combined_weights(n)
            !! Combined weights: product of user weights and robust weights (weights(i) * robust_weights(i))
        real(real64), intent(out) :: tmp_residuals(n)
            !! Residuals (y - fitted_values), used to compute bisquare robust weights
        integer(int32), intent(out) :: tmp_permutation_indices(n)
            !! Permutation indices array (from NetLib bisquare weight computation)

        real(real64), intent(out) :: fitted_values(n)
            !! Fitted (smoothed) values of y at the evaluation points
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: neighborhood_size
            !! Size of the local neighborhood used at the current span (points per local fit)

        integer(int32) :: iter, i, predictor_dim
        integer(int32) :: actual_n_iters
        logical :: actual_compute_influence, actual_save_factorization

        call set_ok(ierr)

        M_DEFAULT_VAL(n_iters, actual_n_iters, CM_DEFAULT_LOESS_ITERS)
        M_DEFAULT_VAL(compute_influence, actual_compute_influence, .false.)
        M_DEFAULT_VAL(save_factorization, actual_save_factorization, .false.)

        ! `span` is a fraction of `n` points included in each local neighborhood; neighborhood_size is that
        ! neighborhood size. Require at least degree+3 points per neighborhood so the local
        ! polynomial fit (degree+1 coefficients) is over-determined rather than exactly/under-determined.
        neighborhood_size = max(2_int32, int(ceiling(span*real(n, real64))))
        call validate_in_range_int(neighborhood_size, ierr, min=degree + 3_int32)
        if (is_err(ierr)) return

        ! Initialize robust weights to 1 (no reweighting on first iteration)
        tmp_robust_weights = 1.0_real64
        predictor_dim = 1_int32

        ! Perform robust iterative refinement
        do iter = 1, actual_n_iters
            do concurrent (i = 1:n) shared(tmp_combined_weights, weights, tmp_robust_weights)
                tmp_combined_weights(i) = weights(i)*tmp_robust_weights(i)
            end do

            ! Perform LOESS fitting for this robust iteration.
            ! `106` is netlib's packed `iv(19)` model-selection code, see loess_fit_plain for details.
            ! the netlib routines take these mandatorily, so hand them the resolved values --
            ! passing the optional dummies straight through crashes whenever the caller omits them
            call loess_decomposition(106_int32, tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace, predictor_dim, n, span, degree, max_neighborhood_size, actual_save_factorization)
            call loess_fitting(x, y, tmp_combined_weights, tmp_hat_diag, actual_compute_influence, tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace)
            call loess_evaluation(tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace, n, eval_points, fitted_values)

            ! Compute residuals for robust reweighting in next iteration
            do concurrent (i = 1:n) shared(tmp_residuals, y, fitted_values)
                tmp_residuals(i) = y(i) - fitted_values(i)
            end do

            ! Compute new robust weights using bisquare reweighting
            call loess_robust_weights(tmp_residuals, n, tmp_robust_weights, tmp_permutation_indices)
        end do
    end subroutine loess_fit_robust_kernel

end module tox_loess_kernel
