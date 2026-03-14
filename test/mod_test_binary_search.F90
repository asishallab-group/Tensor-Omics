#include "macros.h"

!> Unit test suite for binary_search routine.
module mod_test_binary_search
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_utils
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

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(2)

        all_tests(1) = test_case("test_binary_search_insertion", test_binary_search_insertion)
        all_tests(2) = test_case("test_binary_search_insertion", test_binary_search)
    end function get_all_tests

    !> Run all binary_search tests.
    subroutine run_all_tests_binary_search
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All binary_search tests passed successfully."
    end subroutine run_all_tests_binary_search

    !> Run specific binary_search tests by name.
    subroutine run_named_tests_binary_search(test_names)
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
    end subroutine run_named_tests_binary_search

    subroutine test_binary_search_insertion()
        real(real64), dimension(7) :: arr
        integer(int32), dimension(7) :: perm
        integer(int32) :: idx

        ! Sorted values with a duplicate and a NaN at the end
        arr  = [1.0_real64, 3.0_real64, 3.0_real64, 5.0_real64, 7.0_real64, 9.0_real64, M_NAN]
        perm = [1,2,3,4,5,6,7]   ! already sorted, NaN last

        ! Insert before first element
        idx = binary_search_insertion(arr, perm, 0.5_real64)
        call assert_equal_int(idx, 1_int32, "test_binary_search_insertion: Insert before first element")

        ! Insert equal to first element
        idx = binary_search_insertion(arr, perm, 1.0_real64)
        call assert_equal_int(idx, 1_int32, "test_binary_search_insertion: Insert equal to first element")

        ! Insert equal to duplicate value (3.0)
        idx = binary_search_insertion(arr, perm, 3.0_real64)
        call assert_equal_int(idx, 2_int32, "test_binary_search_insertion: Insert equal to duplicate value")

        ! Insert between duplicates
        idx = binary_search_insertion(arr, perm, 4.0_real64)
        call assert_equal_int(idx, 4_int32, "test_binary_search_insertion: Insert between duplicates")

        ! Insert before NaN (should treat NaN as largest)
        idx = binary_search_insertion(arr, perm, 10.0_real64)
        call assert_equal_int(idx, 7_int32, "test_binary_search_insertion: Insert before NaN")

        ! Insert after all real values (NaN is last)
        idx = binary_search_insertion(arr, perm, huge(1.0_real64))
        call assert_equal_int(idx, 7_int32, "test_binary_search_insertion: Insert before NaN with huge value")

        ! Insert using restricted range (lower > upper → should return lower)
        idx = binary_search_insertion(arr, perm, 4.0_real64, lower_idx=5, upper_idx=4)
        call assert_equal_int(idx, 5_int32, "test_binary_search_insertion: Insert with empty search range")

        ! Insert using restricted range inside duplicates
        idx = binary_search_insertion(arr, perm, 3.5_real64, lower_idx=2, upper_idx=3)
        call assert_equal_int(idx, 4_int32, "test_binary_search_insertion: Insert inside restricted duplicate range")
    end subroutine test_binary_search_insertion

    subroutine test_binary_search()
        real(real64), dimension(7) :: arr
        integer(int32), dimension(7) :: perm
        integer(int32) :: idx

        ! Unsorted array with duplicates and NaN
        arr  = [10.0_real64, 5.0_real64, 20.0_real64, 15.0_real64, 5.0_real64, 30.0_real64, M_NAN]
        perm = [2,5,1,4,3,6,7]
        ! arr(perm) = [5,5,10,15,20,30,NaN]

        ! Find first duplicate value
        idx = binary_search(arr, perm, 5.0_real64)
        call assert_equal_int(idx, 1_int32, "test_binary_search: Find first occurrence of duplicate")

        ! Find middle value
        idx = binary_search(arr, perm, 15.0_real64)
        call assert_equal_int(idx, 4_int32, "test_binary_search: Find middle value")

        ! Find last real value
        idx = binary_search(arr, perm, 30.0_real64)
        call assert_equal_int(idx, 6_int32, "test_binary_search: Find last real value")

        ! Search for NaN (should not be found)
        idx = binary_search(arr, perm, M_NAN)
        call assert_equal_int(idx, size(perm, kind=int32), "test_binary_search: NaN should not be found")

        ! Search for value smaller than all elements
        idx = binary_search(arr, perm, -100.0_real64)
        call assert_equal_int(idx, -1_int32, "test_binary_search: Value smaller than all elements")

        ! Search for value larger than all elements
        idx = binary_search(arr, perm, 100.0_real64)
        call assert_equal_int(idx, -1_int32, "test_binary_search: Value larger than all elements")

        ! Search for value between duplicates
        idx = binary_search(arr, perm, 7.0_real64)
        call assert_equal_int(idx, -1_int32, "test_binary_search: Value between duplicates but not present")

        ! Search in array with only NaN
        arr  = M_NAN
        idx = binary_search(arr, perm, 1.0_real64)
        call assert_equal_int(idx, -1_int32, "test_binary_search: Search non-NaN in all-NaN array")
        idx = binary_search(arr, perm, M_NAN)
        call assert_equal_int(idx, size(perm, kind=int32), "test_binary_search: Search non-NaN in all-NaN array")

    end subroutine test_binary_search
end module mod_test_binary_search
