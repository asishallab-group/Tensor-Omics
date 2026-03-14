#include "macros.h"

!> Unit test suite for heaps.
module mod_test_heaps
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use f42_heaps
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 0.0_real64

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(3)

        all_tests(1) = test_case("test_top_k_heap", test_top_k_heap)
        all_tests(2) = test_case("test_min_heap", test_min_heap)
        all_tests(3) = test_case("test_max_heap", test_max_heap)
    end function get_all_tests

    !> Run all heaps tests.
    subroutine run_all_tests_heaps
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All heaps tests passed successfully."
    end subroutine run_all_tests_heaps

    !> Run specific heaps tests by name.
    subroutine run_named_tests_heaps(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_heaps

    subroutine test_top_k_heap()
        integer(int32), parameter :: heap_size = 5
        real(real64), dimension(heap_size) :: heap, expected_heap

        ! Initialize
        call init_top_k_heap(heap, heap_size)
        expected_heap = [M_NEG_INF, M_NEG_INF, M_NEG_INF, M_NEG_INF, M_NEG_INF]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Init")

        ! Push values
        call top_k_heap_push(heap, heap_size, M_NAN)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Push NaN on fresh initializes heap, should be ignored")

        call top_k_heap_push(heap, heap_size, 5.0_real64)
        expected_heap = [M_NEG_INF, M_NEG_INF, M_NEG_INF, 5.0_real64, M_NEG_INF]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: First push")

        call top_k_heap_push(heap, heap_size, 7.0_real64)
        expected_heap = [M_NEG_INF, M_NEG_INF, M_NEG_INF, 5.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Second push")

        call top_k_heap_push(heap, heap_size, 10.0_real64)
        expected_heap = [M_NEG_INF, 5.0_real64, M_NEG_INF, 10.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Third push")

        call top_k_heap_push(heap, heap_size, 6.0_real64)
        expected_heap = [M_NEG_INF, 5.0_real64, 6.0_real64, 10.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Fourth push")

        call top_k_heap_push(heap, heap_size, 11.0_real64)
        expected_heap = [5.0_real64, 7.0_real64, 6.0_real64, 10.0_real64, 11.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Fifth push")

        call top_k_heap_push(heap, heap_size, 4.0_real64)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Push worse, should be ignored")

        call top_k_heap_push(heap, heap_size, 12.0_real64)
        expected_heap = [6.0_real64, 7.0_real64, 12.0_real64, 10.0_real64, 11.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Push on full heap")

        ! Push NaN (should be ignored)
        call top_k_heap_push(heap, heap_size, M_NAN)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_top_k_heap: Push NaN on full heap, should be ignored")
    end subroutine test_top_k_heap

    subroutine test_min_heap()
        integer(int32), parameter :: heap_size = 5
        real(real64), dimension(heap_size) :: heap, expected_heap
        integer(int32) :: n, ierr
        real(real64) :: val

        ! Initialize
        n = 0_int32
        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 3_int32), "test_min_heap: Expected error for popping from empty heap (n=0)")
        n = 1_int32
        call minheap_pop(heap, 0_int32, n, val, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 2_int32), "test_min_heap: Expected error for popping from empty heap (heap_size=0)")

        ! Push values
        n = 0_int32
        heap = M_NAN

        call minheap_push(heap, heap_size, n, M_NAN)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Push NaN on empty heap, should be ignored")
        call assert_equal_int(n, 0_int32, "test_min_heap: Push NaN, expected n unchanged")

        call minheap_push(heap, heap_size, n, 12.0_real64)
        expected_heap = [12.0_real64, M_NAN, M_NAN, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: First push")
        call assert_equal_int(n, 1_int32, "test_min_heap: First push, expected n=1")

        call minheap_push(heap, heap_size, n, 7.0_real64)
        expected_heap = [7.0_real64, 12.0_real64, M_NAN, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Second push")
        call assert_equal_int(n, 2_int32, "test_min_heap: Second push, expected n=2")

        call minheap_push(heap, heap_size, n, 8.0_real64)
        expected_heap = [7.0_real64, 12.0_real64, 8.0_real64, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Third push")
        call assert_equal_int(n, 3_int32, "test_min_heap: Third push, expected n=3")

        call minheap_push(heap, heap_size, n, 6.0_real64)
        expected_heap = [6.0_real64, 7.0_real64, 8.0_real64, 12.0_real64, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fourth push")
        call assert_equal_int(n, 4_int32, "test_min_heap: Fourth push, expected n=4")

        call minheap_push(heap, heap_size, n, 5.0_real64)
        expected_heap = [5.0_real64, 6.0_real64, 8.0_real64, 12.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fifth push")
        call assert_equal_int(n, 5_int32, "test_min_heap: Fifth push, expected n=5")

        call minheap_push(heap, heap_size, n, 5.0_real64)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Push on full heap, should be ignored")
        call assert_equal_int(n, 5_int32, "test_min_heap: Heap already full, expected n unchanged")

        ! Pop values
        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: First pop, ierr")
        call assert_equal_int(n, 4_int32, "test_min_heap: First pop, n")
        call assert_equal_real(val, 5.0_real64, TOL, "test_min_heap: First pop, value")
        expected_heap = [6.0_real64, 7.0_real64, 8.0_real64, 12.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: First pop, expected heap")

        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Second pop, ierr")
        call assert_equal_int(n, 3_int32, "test_min_heap: Second pop, n")
        call assert_equal_real(val, 6.0_real64, TOL, "test_min_heap: Second pop, value")
        expected_heap = [7.0_real64, 12.0_real64, 8.0_real64, 12.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Second pop, expected heap")

        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Third pop, ierr")
        call assert_equal_int(n, 2_int32, "test_min_heap: Third pop, n")
        call assert_equal_real(val, 7.0_real64, TOL, "test_min_heap: Third pop, value")
        expected_heap = [8.0_real64, 12.0_real64, 8.0_real64, 12.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Third pop, expected heap")

        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Fourth pop, ierr")
        call assert_equal_int(n, 1_int32, "test_min_heap: Fourth pop, n")
        call assert_equal_real(val, 8.0_real64, TOL, "test_min_heap: Fourth pop, value")
        expected_heap = [12.0_real64, 12.0_real64, 8.0_real64, 12.0_real64, 7.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fourth pop, expected heap")

        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Fifth pop, ierr")
        call assert_equal_int(n, 0_int32, "test_min_heap: Fifth pop, n")
        call assert_equal_real(val, 12.0_real64, TOL, "test_min_heap: Fifth pop, value")
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fifth pop, expected heap")
    end subroutine test_min_heap

    subroutine test_max_heap()
        integer(int32), parameter :: heap_size = 5
        real(real64), dimension(heap_size) :: heap, expected_heap
        integer(int32) :: n, ierr
        real(real64) :: val

        ! Initialize
        n = 0_int32
        call minheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 3_int32), "test_max_heap: Expected error for popping from empty heap (n=0)")
        n = 1_int32
        call minheap_pop(heap, 0_int32, n, val, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 2_int32), "test_max_heap: Expected error for popping from empty heap (heap_size=0)")

        ! Push values
        n = 0_int32
        heap = M_NAN

        call maxheap_push(heap, heap_size, n, M_NAN)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Push NaN on empty heap, should be ignored")
        call assert_equal_int(n, 0_int32, "test_max_heap: Push NaN, expected n unchanged")

        call maxheap_push(heap, heap_size, n, 5.0_real64)
        expected_heap = [5.0_real64, M_NAN, M_NAN, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: First push")
        call assert_equal_int(n, 1_int32, "test_max_heap: First push, expected n=1")

        call maxheap_push(heap, heap_size, n, 12.0_real64)
        expected_heap = [12.0_real64, 5.0_real64, M_NAN, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Second push")
        call assert_equal_int(n, 2_int32, "test_max_heap: Second push, expected n=2")

        call maxheap_push(heap, heap_size, n, 8.0_real64)
        expected_heap = [12.0_real64, 5.0_real64, 8.0_real64, M_NAN, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Third push")
        call assert_equal_int(n, 3_int32, "test_max_heap: Third push, expected n=3")

        call maxheap_push(heap, heap_size, n, 12.01_real64)
        expected_heap = [12.01_real64, 12.0_real64, 8.0_real64, 5.0_real64, M_NAN]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Fourth push")
        call assert_equal_int(n, 4_int32, "test_max_heap: Fourth push, expected n=4")

        call maxheap_push(heap, heap_size, n, huge(1.0_real64))
        expected_heap = [huge(1.0_real64), 12.01_real64, 8.0_real64, 5.0_real64, 12.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Fifth push")
        call assert_equal_int(n, 5_int32, "test_max_heap: Fifth push, expected n=5")

        call maxheap_push(heap, heap_size, n, 5.0_real64)
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_max_heap: Push on full heap, should be ignored")
        call assert_equal_int(n, 5_int32, "test_max_heap: Heap already full, expected n unchanged")


        ! Pop values
        call maxheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: First pop, ierr")
        call assert_equal_int(n, 4_int32, "test_min_heap: First pop, n")
        call assert_equal_real(val, huge(1.0_real64), TOL, "test_min_heap: First pop, value")
        expected_heap = [12.01_real64, 12.0_real64, 8.0_real64, 5.0_real64, 12.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: First pop, expected heap")

        call maxheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Second pop, ierr")
        call assert_equal_int(n, 3_int32, "test_min_heap: Second pop, n")
        call assert_equal_real(val, 12.01_real64, TOL, "test_min_heap: Second pop, value")
        expected_heap = [12.0_real64, 5.0_real64, 8.0_real64, 5.0_real64, 12.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Second pop, expected heap")

        call maxheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Third pop, ierr")
        call assert_equal_int(n, 2_int32, "test_min_heap: Third pop, n")
        call assert_equal_real(val, 12.0_real64, TOL, "test_min_heap: Third pop, value")
        expected_heap = [8.0_real64, 5.0_real64, 8.0_real64, 5.0_real64, 12.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Third pop, expected heap")

        call maxheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Fourth pop, ierr")
        call assert_equal_int(n, 1_int32, "test_min_heap: Fourth pop, n")
        call assert_equal_real(val, 8.0_real64, TOL, "test_min_heap: Fourth pop, value")
        expected_heap = [5.0_real64, 5.0_real64, 8.0_real64, 5.0_real64, 12.0_real64]
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fourth pop, expected heap")

        call maxheap_pop(heap, heap_size, n, val, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_min_heap: Fifth pop, ierr")
        call assert_equal_int(n, 0_int32, "test_min_heap: Fifth pop, n")
        call assert_equal_real(val, 5.0_real64, TOL, "test_min_heap: Fifth pop, value")
        call assert_equal_array_real(heap, expected_heap, heap_size, TOL, "test_min_heap: Fifth pop, expected heap")
    end subroutine test_max_heap
end module mod_test_heaps
