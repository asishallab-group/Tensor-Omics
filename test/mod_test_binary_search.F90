#include <src/macros.h>

!> Unit test suite for binary_search routine.
module mod_test_binary_search
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_utils_impl
    use tox_errors
    use test_suite, only: test_case
    implicit none

contains

    !> Get array of all available tests.
    function get_all_tests_binary_search() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(3))

        all_tests(1) = test_case("test_binary_search_insertion", test_binary_search_insertion)
        all_tests(2) = test_case("test_binary_search", test_binary_search)
        all_tests(3) = test_case("test_block_containing_index", test_block_containing_index)
    end function get_all_tests_binary_search

    subroutine test_binary_search_insertion()
        real(real64), dimension(7) :: arr
        integer(int32), dimension(7) :: perm
        integer(int32) :: idx

        ! Sorted values with a duplicate and a NaN at the end
        arr = [1.0_real64, 3.0_real64, 3.0_real64, 5.0_real64, 7.0_real64, 9.0_real64, M_NAN]
        perm = [1, 2, 3, 4, 5, 6, 7]   ! already sorted, NaN last

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
        arr = [10.0_real64, 5.0_real64, 20.0_real64, 15.0_real64, 5.0_real64, 30.0_real64, M_NAN]
        perm = [2, 5, 1, 4, 3, 6, 7]
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
        arr = M_NAN
        idx = binary_search(arr, perm, 1.0_real64)
        call assert_equal_int(idx, -1_int32, "test_binary_search: Search non-NaN in all-NaN array")
        idx = binary_search(arr, perm, M_NAN)
        call assert_equal_int(idx, size(perm, kind=int32), "test_binary_search: Search non-NaN in all-NaN array")

    end subroutine test_binary_search

    !> block_containing_index on block ends `[3, 3, 4, 8, 8]` (block 2 and block 5 empty): every
    !| index 1..8 lands in the first block whose end reaches it, so an empty block is never chosen
    !| (index 4 skips block 2 for block 3, index 8 stays in block 4 rather than the empty block
    !| 5). The first and last index exercise both ends of the search; a single block holds
    !| everything.
    subroutine test_block_containing_index()
        character(len=*), parameter :: name = "test_block_containing_index: "
        integer(int32), parameter :: ends(5) = [3_int32, 3_int32, 4_int32, 8_int32, 8_int32]
        integer(int32), parameter :: expected(8) = [1_int32, 1_int32, 1_int32, 3_int32, 4_int32, 4_int32, 4_int32, &
                                                    4_int32]
        integer(int32) :: idx, got(8)

        do idx = 1, 8
            got(idx) = block_containing_index(ends, 5_int32, idx)
        end do
        call assert_equal_array_int(got, expected, 8_int32, name//"five blocks with two empty ones")

        call assert_equal_int(block_containing_index([6_int32], 1_int32, 1_int32), 1_int32, &
                              name//"one block, first index")
        call assert_equal_int(block_containing_index([6_int32], 1_int32, 6_int32), 1_int32, &
                              name//"one block, last index")
        call assert_equal_int(block_containing_index([1_int32, 2_int32], 2_int32, 2_int32), 2_int32, &
                              name//"two blocks of one index, second")
    end subroutine test_block_containing_index
end module mod_test_binary_search
