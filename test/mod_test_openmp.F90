!> Unit test suite for f42_openmp: what a build reports about OpenMP.
!> The same tests run in both builds; the `!$` lines are compiled only with OpenMP.
!> What the build was MEANT to be comes from TOX_EXPECT_OPENMP, which test_runner.sh sets from its
!> `--openmp` option (1 with it, 0 without). Run by hand without the variable, the test expects
!> what it was compiled as -- which a build that lost or gained OpenMP by mistake agrees with.

module mod_test_openmp
  use asserts
  use f42_openmp, only: max_omp_threads
  use, intrinsic :: iso_fortran_env, only: int32
  use test_suite
  implicit none
  public

contains

  !> Get array of all available OpenMP tests.
  function get_all_tests_openmp() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(1))

    all_tests(1) = test_case("test_max_omp_threads_follows_the_runtime", test_max_omp_threads_follows_the_runtime)
  end function get_all_tests_openmp

  !> A serial build answers 1 whatever the runtime is told; an OpenMP build answers what it is told.
  !> With TOX_EXPECT_OPENMP=1 the answer must be the requested count, with 0 it must be 1.
  subroutine test_max_omp_threads_follows_the_runtime()
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    integer(int32), parameter :: n_threads_requested = 3_int32
    character(len=8) :: expect_openmp
    integer :: length, status
    integer(int32) :: expected, n_threads_before

    ! unset (a run by hand): what this file was compiled as
    expected = 1_int32
    n_threads_before = 1_int32
    !$ n_threads_before = int(omp_get_max_threads(), int32)
    !$ call omp_set_num_threads(n_threads_requested)
    !$ expected = n_threads_requested

    ! set (by the test runner): what the build was meant to be, whatever it was compiled as
    call get_environment_variable("TOX_EXPECT_OPENMP", expect_openmp, length, status)
    if (length > 0) then
      if (status == 0 .and. expect_openmp == "1") then
        expected = n_threads_requested
      else if (status == 0 .and. expect_openmp == "0") then
        expected = 1_int32
      else
        call assert_true(.false., "TOX_EXPECT_OPENMP is 1, 0 or unset, not '" // trim(expect_openmp) // "'")
      end if
    end if

    call assert_equal_int(max_omp_threads(), expected, "max_omp_threads()")
    !$ call omp_set_num_threads(n_threads_before)
  end subroutine test_max_omp_threads_follows_the_runtime

end module mod_test_openmp
