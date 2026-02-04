#include "macros.h"

module tox_loess
  use safeguard
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
  
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
  interface
    ! ============================================================
    ! Subroutine: lowesd
    ! ============================================================
    !> Perform LOESS decomposition.
    !| This subroutine computes the decomposition of the LOESS model.
    !| It initializes the integer and real workspace arrays based on the input parameters.
    subroutine lowesd(mode, iv, liv, lv, wv, dim, n, span, degree, nvmax, setlf)
      use, intrinsic :: iso_fortran_env, only: real64, int32
      integer(int32), intent(in) :: mode
      !| Mode of operation
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      integer(int32), intent(in) :: dim
      !| Dimensionality of the data
      integer(int32), intent(in) :: n
      !| Number of data points
      integer(int32), intent(in) :: degree
      !| Degree of the LOESS polynomial
      integer(int32), intent(in) :: nvmax
      !| Maximum neighborhood size
      real(real64), intent(in) :: span
      !| Smoothing parameter for LOESS
      logical, intent(in) :: setlf
      !| Save matrix factorization flag
      integer(int32), intent(inout) :: iv(liv)
      !| Integer workspace array
      real(real64), intent(inout) :: wv(lv)
      !| Real workspace array
    end subroutine lowesd

    ! ============================================================
    ! Subroutine: lowesb
    ! ============================================================
    !> Perform LOESS fitting.
    !| This subroutine fits the LOESS model to the data.
    !| It uses the input predictor and response variables to compute the diagonal elements of the hat matrix.
    subroutine lowesb(x, y, w, diagl, infl, iv, liv, lv, wv)
      import :: real64, int32
      real(real64), intent(in) :: x(*)
      !| Predictor variable array
      real(real64), intent(in) :: y(*)
      !| Response variable array
      real(real64), intent(in) :: w(*)
      !| Weight array for data points
      real(real64), intent(out) :: diagl(*)
      !| Diagonal elements of the hat matrix
      logical, intent(in) :: infl
      !| Influence calculation flag
      integer(int32), intent(in) :: iv(*)
      !| Integer workspace array
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      real(real64), intent(inout) :: wv(*)
      !| Real workspace array
    end subroutine lowesb

    ! ============================================================
    ! Subroutine: lowese
    ! ============================================================
    !> Perform LOESS evaluation.
    !| This subroutine evaluates the LOESS model.
    !| It computes the smoothed response variable array based on the input predictor variables.
    subroutine lowese(iv, liv, lv, wv, n, z, s)
      import :: real64, int32
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      integer(int32), intent(in) :: n
      !| Number of data points
      integer(int32), intent(inout) :: iv(*)
      !| Integer workspace array
      real(real64), intent(inout) :: wv(*)
      !| Real workspace array
      real(real64), intent(in) :: z(*)
      !| Predictor variable array
      real(real64), intent(out) :: s(*)
      !| Smoothed response variable array
    end subroutine lowese

    ! ============================================================
    ! Subroutine: lowesw
    ! ============================================================
    !> Compute robust weights for LOESS.
    !| This subroutine computes the robust weights for the LOESS model.
    !| It uses the residuals to update the weights for robust fitting.
    subroutine lowesw(res, n, rw, pi)
      use, intrinsic :: iso_fortran_env, only: real64, int32
      real(real64), intent(in) :: res(*)
      !| Residuals array
      integer(int32), intent(in) :: n
      !| Number of data points
      real(real64), intent(out) :: rw(*)
      !| Robust weights array
      integer(int32), intent(out) :: pi(*)
      !| Permutation indices array
    end subroutine lowesw
  end interface

contains

  ! ============================================================
  ! Recommend workspace sizes based on Netlib exact formulas
  ! ============================================================
  !> Recommend workspace sizes based on Netlib exact formulas.
  !| Computes the required sizes for integer and real workspace arrays.
  !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
  subroutine tox_loess_required_workspace(d, nvmax, liv, lv, setlf)
    !! Recommend sizes for integer pool iv(:) and real pool v(:).
    !! nvmax is the maximum neighborhood size (usually n_train).
    !! setlf indicates if matrix factorizations are saved.
    integer(int32), intent(in) :: d
    !| Dimensionality of the data
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    logical, intent(in) :: setlf
    !| Save matrix factorization flag
    integer(int32), intent(out) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(out) :: lv
    !| Length of the real workspace array

    ! Adjusted formulas to ensure sufficient workspace
    liv = 100 + (2**d + 15) * nvmax 
    if (setlf) liv = liv + nvmax
    lv = 100 + (10 * d + 20) * nvmax
    if (setlf) lv = lv + (d + 1) * nvmax
    liv = max(10000_int32, liv)
    lv  = max(100000_int32, lv)
    ! liv  = max(50000, 500 * nvmax)
    ! lv   = max(500000, 5000 * nvmax)
  end subroutine tox_loess_required_workspace

  ! ============================================================
  ! Plain LOESS fitting
  ! ============================================================
  !> Perform plain LOESS fitting.
  !| Fits a LOESS model to the data using the specified smoothing parameter.
  !| Outputs the smoothed response variable array.
  subroutine loess_fit_plain(n, x, y, w, z, span, degree, nvmax, infl, setlf, iv, liv, wv, lv, diagl, yhat, ierr)
    integer(int32), intent(in) :: n
    !| Total number of data points
    integer(int32), intent(in) :: degree
    !| Degree of the LOESS polynomial
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    integer(int32), intent(in) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(in) :: lv
    !| Length of the real workspace array

    real(real64), intent(in) :: x(n)
    !| Predictor variable array
    real(real64), intent(in) :: y(n)
    !| Response variable array
    real(real64), intent(in) :: w(n)
    !| Weight array for data points
    real(real64), intent(in) :: z(n,1)
    !| Additional predictor variable array
    real(real64), intent(in) :: span
    !| Smoothing parameter for LOESS

    logical, intent(in) :: infl
    !| Influence calculation flag
    logical, intent(in) :: setlf
    !| Save matrix factorization flag

    integer(int32), intent(inout) :: iv(liv)
    !| Integer workspace array
    real(real64), intent(inout) :: wv(lv)
    !| Real workspace array
    real(real64), intent(inout) :: diagl(n)
    !| Diagonal elements of the hat matrix

    real(real64), intent(out) :: yhat(n)
    !| Smoothed response variable array
    integer(int32), intent(out) :: ierr
    !! Error code

    ! Initialize error code
    call set_ok(ierr)

    ! Validate inputs
    call validate_dimension_size(n, ierr)
    if (is_err(ierr)) then
      return
    end if

    if (n == 1) then
      yhat(1) = y(1)
      call set_err(ierr, ERR_INVALID_INPUT)
      return 
    end if

    call validate_in_range_real(span, ierr, min=0.0_real64)
    if (is_err(ierr)) then
      return
    end if

    ! Validate workspace sizes
    if (liv < 10000 .or. lv < 100000) then
      call set_err(ierr, ERR_INVALID_INPUT)
      return
    end if

    call lowesd(106, iv, liv, lv, wv, 1, n, span, degree, nvmax, setlf)
    call lowesb(x, y, w, diagl, infl, iv, liv, lv, wv)
    call lowese(iv, liv, lv, wv, n, z, yhat)

  end subroutine loess_fit_plain

  ! ============================================================
  ! Robust LOESS fitting 
  ! ============================================================
  !> Perform robust LOESS fitting with bisquare reweighting.
  !| Fits a LOESS model to the data using robust iterations to handle outliers.
  !| Outputs the smoothed response variable array.
  subroutine loess_fit_robust(n, x, y, w, z, span, degree, nvmax, infl, setlf, n_iters, iv, liv, wv, lv, diagl, rw, ww, res, pi, yhat, ierr)
    integer(int32), intent(in) :: n
    !| Total number of data points
    integer(int32), intent(in) :: degree
    !| Degree of the LOESS polynomial
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    integer(int32), intent(in) :: n_iters
    !| Number of robust iterations
    integer(int32), intent(in) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(in) :: lv
    !| Length of the real workspace array

    real(real64), intent(in) :: x(n)
    !| Predictor variable array
    real(real64), intent(in) :: y(n)
    !| Response variable array
    real(real64), intent(in) :: w(n)
    !| Weight array for data points
    real(real64), intent(in) :: z(n,1)
    !| Additional predictor variable array
    real(real64), intent(in) :: span
    !| Smoothing parameter for LOESS

    logical, intent(in) :: infl
    !| Influence calculation flag
    logical, intent(in) :: setlf
    !| Save matrix factorization flag

    integer(int32), intent(inout) :: iv(liv)
    !| Integer workspace array
    real(real64), intent(inout) :: wv(lv)
    !| Real workspace array
    real(real64), intent(inout) :: diagl(n)
    !| Diagonal elements of the hat matrix
    real(real64), intent(inout) :: rw(n)
    !| Robust weights array
    real(real64), intent(inout) :: ww(n)
    !| Working weights array
    real(real64), intent(inout) :: res(n)
    !| Residuals array
    integer(int32), intent(inout) :: pi(n)
    !| Permutation indices array

    real(real64), intent(out) :: yhat(n)
    !| Smoothed response variable array
    integer(int32), intent(out) :: ierr
    !! Error code

    integer :: it, i, dim_val

    ! Initialize error code
    call set_ok(ierr)

    ! Validate inputs
    call validate_dimension_size(n, ierr)
    call validate_in_range_real(span, ierr, min=0.0_real64)
    if (is_err(ierr)) return

    ! Validate workspace sizes
    if (liv < 10000 .or. lv < 100000) then
      call set_err(ierr, ERR_INVALID_INPUT)
      return
    end if

    if (n == 1) then
      yhat(1) = y(1)
      call set_err(ierr, ERR_INVALID_INPUT)
      return 
    end if

    rw = 1.0_real64
    dim_val = 1

    do it = 1, n_iters
      do i = 1, n
        ww(i) = w(i) * rw(i)
      end do

      call lowesd(106, iv, liv, lv, wv, dim_val, n, span, degree, nvmax, setlf)
      call lowesb(x, y, ww, diagl, infl, iv, liv, lv, wv)
      call lowese(iv, liv, lv, wv, n, z, yhat)

      do i = 1, n
        res(i) = y(i) - yhat(i)
      end do

      call lowesw(res, n, rw, pi)
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
  subroutine loess_alloc(x, y, span, degree, yhat, mode, n_iters, ierr)
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Input parameters
    real(real64), intent(in) :: x(:)
    !| Predictor variable array
    real(real64), intent(in) :: y(:)
    !| Response variable array
    real(real64), intent(in) :: span
    !| Smoothing parameter for LOESS
    integer(int32), intent(in) :: degree
    !| Degree of the LOESS polynomial
    integer(int32), intent(in) :: mode
    !| Mode of operation: 0 for plain, 1 for robust
    integer(int32), intent(in) :: n_iters
    !| Number of robust iterations (only used when mode = 1)

    ! Output parameters
    real(real64), intent(out) :: yhat(size(y))
    !| Smoothed response variable array
    integer(int32), intent(out) :: ierr
    !| Error code

    ! Local variables
    integer(int32) :: n, liv, lv
    integer(int32), allocatable :: iv(:), pi(:)
    real(real64), allocatable :: wv(:), diagl(:), rw(:), ww(:), res(:), w_init(:), z_mat(:,:)

    ! Initialize variables
    n = size(y)
    call set_ok(ierr)

    ! Compute workspace sizes based on Netlib formulas
    call tox_loess_required_workspace(1_int32, n, liv, lv, .false.)

    ! Allocate workspace arrays
    allocate(iv(liv), wv(lv), diagl(n), w_init(n), z_mat(n, 1), stat=ierr)
    if (ierr /= 0) return

    ! Initialize arrays
    iv = 1_int32
    w_init = 1.0_real64
    z_mat(:, 1) = x

    ! Allocate additional arrays for robust mode
    if (mode == 1) then
      allocate(rw(n), ww(n), res(n), pi(n), stat=ierr)
      if (ierr /= 0) then
        deallocate(iv, wv, diagl, w_init, z_mat)
        return
      end if
    end if

    ! Call the appropriate LOESS fitting subroutine
    if (mode == 0) then
      call loess_fit_plain(n, x, y, w_init, z_mat, span, degree, n, &
                          .false., .false., iv, liv, wv, lv, diagl, yhat, ierr)
    else
      call loess_fit_robust(n, x, y, w_init, z_mat, span, degree, n, &
                            .false., .false., n_iters, iv, liv, wv, lv, &
                            diagl, rw, ww, res, pi, yhat, ierr)
    end if

    ! Deallocate workspace arrays
    deallocate(iv, wv, diagl, w_init, z_mat)
    if (allocated(rw)) deallocate(rw, ww, res, pi)

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
subroutine tox_loess_c(x, y, n, span, degree, yhat, mode, n_iters, ierr) bind(C, name="tox_loess_c")
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
  real(c_double), dimension(n), intent(out), target :: yhat
  !! Smoothed response variable array
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
  M_CHECK_NON_NULL(yhat)

  ! Call the Fortran subroutine
  call loess_alloc(x, y, span, degree, yhat, mode, n_iters, ierr)
end subroutine tox_loess_c

!> Perform plain LOESS fitting.
!| Fits a LOESS model to the data using the specified smoothing parameter.
!| Outputs the smoothed response variable array.
subroutine loess_fit_plain_c(n, x, y, w, z, span, degree, nvmax, infl, setlf, iv, liv, wv, lv, diagl, yhat, ierr) bind(C, name="loess_fit_plain_c")
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
  real(c_double), dimension(n, 1), intent(in), target :: z
  !! Additional predictor variable array
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
  integer(c_int), dimension(*), intent(inout), target :: iv
  !! Integer workspace array
  integer(c_int), intent(in), target :: liv
  !! Length of the integer workspace array
  real(c_double), dimension(*), intent(inout), target :: wv
  !! Real workspace array
  integer(c_int), intent(in), target :: lv
  !! Length of the real workspace array

  ! Output parameters
  real(c_double), dimension(n), intent(out), target :: diagl
  !! Diagonal elements of the hat matrix
  real(c_double), dimension(n), intent(out), target :: yhat
  !! Smoothed response variable array
  integer(c_int), intent(out), target :: ierr
  !! Error code

  ! Null validation macros
  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n)
  M_CHECK_NON_NULL(x)
  M_CHECK_NON_NULL(y)
  M_CHECK_NON_NULL(w)
  M_CHECK_NON_NULL(z)
  M_CHECK_NON_NULL(span)
  M_CHECK_NON_NULL(degree)
  M_CHECK_NON_NULL(nvmax)
  M_CHECK_NON_NULL(infl)
  M_CHECK_NON_NULL(setlf)
  M_CHECK_NON_NULL(iv)
  M_CHECK_NON_NULL(liv)
  M_CHECK_NON_NULL(wv)
  M_CHECK_NON_NULL(lv)
  M_CHECK_NON_NULL(diagl)
  M_CHECK_NON_NULL(yhat)

  ! Call the Fortran subroutine
  call loess_fit_plain(n, x, y, w, z, span, degree, nvmax, infl /= 0, setlf /= 0, iv, liv, wv, lv, diagl, yhat, ierr)
end subroutine loess_fit_plain_c

!> Perform robust LOESS fitting with bisquare reweighting.
!| Fits a LOESS model to the data using robust iterations to handle outliers.
!| Outputs the smoothed response variable array.
subroutine loess_fit_robust_c(n, x, y, w, z, span, degree, nvmax, infl, setlf, n_iters, iv, liv, wv, lv, diagl, rw, ww, res, pi, yhat, ierr) bind(C, name="loess_fit_robust_c")
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
  real(c_double), dimension(n, 1), intent(in), target :: z
  !! Additional predictor variable array
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
  integer(c_int), dimension(*), intent(inout), target :: iv
  !! Integer workspace array
  integer(c_int), intent(in), target :: liv
  !! Length of the integer workspace array
  real(c_double), dimension(*), intent(inout), target :: wv
  !! Real workspace array
  integer(c_int), intent(in), target :: lv
  !! Length of the real workspace array

  ! Output parameters
  real(c_double), dimension(n), intent(out), target :: diagl
  !! Diagonal elements of the hat matrix
  real(c_double), dimension(n), intent(out), target :: rw
  !! Robust weights array
  real(c_double), dimension(n), intent(out), target :: ww
  !! Working weights array
  real(c_double), dimension(n), intent(out), target :: res
  !! Residuals array
  integer(c_int), dimension(n), intent(out), target :: pi
  !! Permutation indices array
  real(c_double), dimension(n), intent(out), target :: yhat
  !! Smoothed response variable array
  integer(c_int), intent(out), target :: ierr
  !! Error code

  ! Null validation macros
  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n)
  M_CHECK_NON_NULL(x)
  M_CHECK_NON_NULL(y)
  M_CHECK_NON_NULL(w)
  M_CHECK_NON_NULL(z)
  M_CHECK_NON_NULL(span)
  M_CHECK_NON_NULL(degree)
  M_CHECK_NON_NULL(nvmax)
  M_CHECK_NON_NULL(infl)
  M_CHECK_NON_NULL(setlf)
  M_CHECK_NON_NULL(n_iters)
  M_CHECK_NON_NULL(iv)
  M_CHECK_NON_NULL(liv)
  M_CHECK_NON_NULL(wv)
  M_CHECK_NON_NULL(lv)
  M_CHECK_NON_NULL(diagl)
  M_CHECK_NON_NULL(rw)
  M_CHECK_NON_NULL(ww)
  M_CHECK_NON_NULL(res)
  M_CHECK_NON_NULL(pi)
  M_CHECK_NON_NULL(yhat)

  ! Call the Fortran subroutine
  call loess_fit_robust(n, x, y, w, z, span, degree, nvmax, infl /= 0, setlf /= 0, n_iters, iv, liv, wv, lv, diagl, rw, ww, res, pi, yhat, ierr)
end subroutine loess_fit_robust_c

!! Wrapper for recommending workspace sizes based on Netlib exact formulas.
!! This subroutine is designed to be called from C code and computes the required sizes
!! for integer and real workspace arrays. These sizes depend on the dimensionality of the data,
!! the maximum neighborhood size, and whether matrix factorizations are saved.
subroutine tox_loess_required_workspace_c(d, nvmax, liv, lv, setlf) bind(C, name="tox_loess_required_workspace_c")
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
  integer(c_int), intent(out) :: liv
  !! Length of the integer workspace array
  integer(c_int), intent(out) :: lv
  !! Length of the real workspace array

  ! Call the Fortran subroutine
  call tox_loess_required_workspace(d, nvmax, liv, lv, setlf /= 0)
end subroutine tox_loess_required_workspace_c