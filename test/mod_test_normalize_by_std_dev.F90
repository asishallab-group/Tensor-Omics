!> Unit test suite for normalize_by_std_dev routine.
module mod_test_normalize_by_std_dev
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite
  use tox_errors
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_normalize_by_std_dev() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(3))
    
    all_tests(1) = test_case("test_std_dev", test_std_dev)
    all_tests(2) = test_case("test_loess_normalization_outlier_correction", test_loess_normalization_outlier_correction)
    all_tests(3) = test_case("test_loess_zero_variance_handling", test_loess_zero_variance_handling)

  end function get_all_tests_normalize_by_std_dev

  subroutine test_std_dev()
    real(real64), dimension(5) :: v
    real(real64) :: s_pop, s_samp

    ! Test vector
    v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

    ! --- Expected values ---
    ! Population variance = 2.0  → std = sqrt(2)
    ! Sample variance     = 2.5  → std = sqrt(2.5)

    ! --- Compute using the user's function ---
    s_pop  = std_dev(v)                     ! default: no Bessel correction
    s_samp = std_dev(v, do_bessel_correction=.true.)

    ! --- Assertions ---
    call assert_equal_real(s_pop,  sqrt(2.0_real64), 1e-12_real64, &
         "test_std_dev: population std_dev incorrect")

    call assert_equal_real(s_samp, sqrt(2.5_real64), 1e-12_real64, &
         "test_std_dev: sample std_dev (Bessel) incorrect")

    ! No NaNs or Infs
    call assert_true(.not. ieee_is_nan(s_pop),  "test_std_dev: population std_dev produced NaN")
    call assert_true(.not. ieee_is_nan(s_samp), "test_std_dev: sample std_dev produced NaN")
    call assert_in_range_real(s_pop,  0.0_real64, huge(1.0_real64), "test_std_dev: population std_dev out of range")
    call assert_in_range_real(s_samp, 0.0_real64, huge(1.0_real64), "test_std_dev: sample std_dev out of range")

    ! Edge case: vector of length 1 → std = 0 for both modes
    block
      real(real64), dimension(1) :: w
      real(real64) :: s1, s1b
      w = [3.14159_real64]

      s1  = std_dev(w)
      s1b = std_dev(w, do_bessel_correction=.true.)

      call assert_equal_real(s1,  0.0_real64, 1e-12_real64, "test_std_dev: std_dev length-1 (population)")
      call assert_equal_real(s1b, 0.0_real64, 1e-12_real64, "test_std_dev: std_dev length-1 (sample)")
    end block

  end subroutine test_std_dev

  !> Main test: Verifies that an SD outlier is corrected by the global curve.
  subroutine test_loess_normalization_outlier_correction()
    integer(int32), parameter :: ng = 20, nt = 10
    real(real64) :: mat(nt, ng), res(nt, ng)
    integer(int32) :: ierr, i
    real(real64) :: span = 0.75_real64
    integer(int32) :: deg = 1

    ! 1. Create data where SD = Mean (Perfect relationship)
    do i = 1, ng
       ! Gene i has mean i, and values fluctuate to give SD approximately i
       mat(:, i) = real(i, real64) 
       mat(1, i) = mat(1, i) + real(i, real64) * 0.5_real64
       mat(2, i) = mat(2, i) - real(i, real64) * 0.5_real64
    end do

    ! 2. Introduce an OUTLIER in gene 10
    ! Assign it a huge variance that breaks the linear trend
    mat(1, 10) = 1000.0_real64 

    call normalize_by_std_dev_alloc(ng, nt, mat, res, span, deg, ierr)

    call assert_equal_int(ierr, ERR_OK, "test_loess_normalization_outlier_correction: LOESS normalization failed")

    call assert_no_nan_real(res, ng*nt, "test_loess_normalization_outlier_correction: NaNs in LOESS result")
  end subroutine test_loess_normalization_outlier_correction

    !> Verifies that genes with SD=0 do not break the routine.
    subroutine test_loess_zero_variance_handling()
        integer(int32), parameter :: ng = 10, nt = 5
        real(real64) :: mat(nt, ng), res(nt, ng)
        integer(int32) :: ierr, i
        
        mat = 1.0_real64 

        do i = 1, 7
          mat(:, i) = 1.0_real64 + real(i, real64) * 0.1_real64
          mat(1, i) = mat(1, i) + 0.5_real64
        end do

        call normalize_by_std_dev_alloc(ng, nt, mat, res, 0.75d0, 1_int32, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_loess_zero_variance_handling: LOESS failed even with valid points")
        
        call assert_equal_real(res(1, 10), 1.0_real64, 1d-12, "test_loess_zero_variance_handling: Zero variance gene altered")
    end subroutine test_loess_zero_variance_handling

end module mod_test_normalize_by_std_dev