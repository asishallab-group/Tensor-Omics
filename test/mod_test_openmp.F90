!> Unit test suite for f42_openmp: what a build reports about OpenMP.
!> The same tests run in both builds; the `!$` lines are compiled only with OpenMP.

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
  subroutine test_max_omp_threads_follows_the_runtime()
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    integer(int32) :: expected, n_threads_before

    expected = 1_int32
    n_threads_before = 1_int32
    !$ n_threads_before = int(omp_get_max_threads(), int32)
    !$ call omp_set_num_threads(3)
    !$ expected = 3_int32
    call assert_equal_int(max_omp_threads(), expected, "max_omp_threads()")
    !$ call omp_set_num_threads(n_threads_before)
  end subroutine test_max_omp_threads_follows_the_runtime

end module mod_test_openmp
