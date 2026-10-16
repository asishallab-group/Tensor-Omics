!> Unit test suite for the f42_sort_impl module.
!| Every assertion helper takes explicit-shape dummies, so a vector subscript like `data(perm)`
!| cannot be passed directly: the compiler has to gather it into an array temporary first, and
!| `-check arg_temp_created` reports every one of them, which buries the test output. The sorted
!| values are therefore gathered into a declared local and that local is what gets asserted on.
module mod_test_sorting
    use f42_utils_impl
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_sorting() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(27))

        all_tests(1) = test_case("test_sort_real", test_sort_real)
        all_tests(2) = test_case("test_sort_integer", test_sort_integer)
        all_tests(3) = test_case("test_sort_character", test_sort_character)
        all_tests(4) = test_case("test_sort_integer_ascending", test_sort_integer_ascending)
        all_tests(5) = test_case("test_sort_real_descending", test_sort_real_descending)
        all_tests(6) = test_case("test_sort_char_random", test_sort_char_random)
        all_tests(7) = test_case("test_sort_sorted_stability", test_sort_sorted_stability)
        all_tests(8) = test_case("test_sort_empty_array", test_sort_empty_array)
        all_tests(9) = test_case("test_sort_large_random", test_sort_large_random)
        ! Heapsort tests (exercise heapsort implementations in f42_sort_impl)
        all_tests(10) = test_case("test_heapsort_real", test_heapsort_real)
        all_tests(11) = test_case("test_heapsort_integer", test_heapsort_integer)
        all_tests(12) = test_case("test_heapsort_character", test_heapsort_character)
        all_tests(13) = test_case("test_heapsort_integer_ascending", test_heapsort_integer_ascending)
        all_tests(14) = test_case("test_heapsort_real_descending", test_heapsort_real_descending)
        all_tests(15) = test_case("test_heapsort_char_random", test_heapsort_char_random)
        all_tests(16) = test_case("test_heapsort_sorted_stability", test_heapsort_sorted_stability)
        all_tests(17) = test_case("test_heapsort_empty_array", test_heapsort_empty_array)
        all_tests(18) = test_case("test_heapsort_large_random", test_heapsort_large_random)
        all_tests(19) = test_case("test_sort_duplicates", test_sort_duplicates)
        all_tests(20) = test_case("test_sort_negatives", test_sort_negatives)
        all_tests(21) = test_case("test_sort_nan", test_sort_nan)
        all_tests(22) = test_case("test_sort_duplicates_real", test_sort_duplicates_real)
        all_tests(23) = test_case("test_sort_negatives_real", test_sort_negatives_real)
        ! Value heaps tracking the k largest / k smallest values pushed
        all_tests(24) = test_case("test_top_k_heap_push", test_top_k_heap_push)
        all_tests(25) = test_case("test_bottom_k_heap_push", test_bottom_k_heap_push)
        all_tests(26) = test_case("test_k_heap_push_size_zero_and_one", test_k_heap_push_size_zero_and_one)
        all_tests(27) = test_case("test_heap_sift_down", test_heap_sift_down)
    end function get_all_tests_sorting

    !> Test sorting of a real array using permutation vector.
    subroutine test_sort_real()
        real(real64), dimension(5) :: data = [3.0d0, 1.0d0, 5.0d0, 2.0d0, 4.0d0]
        integer(int32), dimension(5) :: perm, expected = [2, 4, 1, 5, 3]
        integer(int32) :: stack_left(20), stack_right(20), i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)

        call assert_equal_array_int(perm, expected, 5, "test_sort_real: permutation mismatch")
    end subroutine test_sort_real

    !> Test sorting of an integer array using permutation vector.
    subroutine test_sort_integer()
        integer(int32), dimension(4) :: data = [10, 3, 7, 1]
        integer(int32), dimension(4) :: perm, expected = [4, 2, 3, 1]
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 4)]
        call sort_array(data, perm, stack_left, stack_right)

        call assert_equal_array_int(perm, expected, 4, "test_sort_integer: permutation mismatch")
    end subroutine test_sort_integer

    !> Test sorting of a character array using permutation vector.
    subroutine test_sort_character()
        character(len=6), dimension(3) :: data = ['delta ', 'alpha ', 'beta  ']
        integer(int32), dimension(3) :: perm, expected = [2, 3, 1]
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 3)]
        call sort_array(data, perm, stack_left, stack_right)

        call assert_equal_array_int(perm, expected, 3, "test_sort_character: permutation mismatch")
    end subroutine test_sort_character

    !> Test that sorted values are in ascending order for integers.
    subroutine test_sort_integer_ascending()
        integer(int32), dimension(5) :: data = [5, 2, 9, 1, 6]
        integer(int32), dimension(5) :: perm, expected_sorted = [1, 2, 5, 6, 9]
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)

        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 5, "test_sort_integer_ascending: sorted values mismatch")
    end subroutine test_sort_integer_ascending

    !> Test that sorted values are in ascending order for reals.
    subroutine test_sort_real_descending()
        real(real64), dimension(5) :: data = [3.5d0, 2.2d0, 8.8d0, 1.1d0, 7.7d0]
        real(real64), dimension(5) :: expected_sorted = [1.1d0, 2.2d0, 3.5d0, 7.7d0, 8.8d0]
        real(real64), dimension(5) :: sorted_data
        integer(int32), dimension(5) :: perm
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)

        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 5, 1d-12, &
                                     "test_sort_real_descending: sorted values mismatch")
    end subroutine test_sort_real_descending

    !> Test sorting of character array with lexicographic ordering.
    subroutine test_sort_char_random()
        character(len=8), dimension(5) :: data = ['dog     ', 'apple   ', 'zebra   ', 'cat     ', 'bird    ']
        character(len=8), dimension(5) :: expected = ['apple   ', 'bird    ', 'cat     ', 'dog     ', 'zebra   ']
        integer(int32), dimension(5) :: perm
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)

        call assert_true(all(data(perm) == expected), "test_sort_char_random: sorted values mismatch")
    end subroutine test_sort_char_random

    !> Test sorting stability with already sorted array.
    subroutine test_sort_sorted_stability()
        integer(int32), dimension(5) :: data = [1, 2, 3, 4, 5]
        integer(int32), dimension(5) :: expected = [1, 2, 3, 4, 5]
        integer(int32), dimension(5) :: perm
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)

        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected, 5, "test_sort_sorted_stability: sorted values mismatch")
    end subroutine test_sort_sorted_stability

    !> Test sorting with empty array (edge case).
    subroutine test_sort_empty_array()
        integer(int32), dimension(0) :: data
        integer(int32), dimension(0) :: perm
        integer(int32) :: stack_left(1), stack_right(1)

        ! This should not crash - just testing robustness
        call sort_array(data, perm, stack_left, stack_right)

        call assert_true(.true., "test_sort_empty_array: empty array handling")
    end subroutine test_sort_empty_array

    !> Test sorting with large random array for performance and correctness.
    subroutine test_sort_large_random()
        real(real64), allocatable :: rdata(:)
        integer(int32), allocatable :: data(:), perm(:), sorted(:), sorted_data(:)
        integer(int32) :: n
        integer(int32), allocatable :: stack_left(:), stack_right(:), dummy_perm(:)
        integer(int32) :: i
        integer(int32) :: n_seed
        integer(int32), allocatable :: seed_array(:)
        n = 1000
        allocate (rdata(n), data(n), perm(n), sorted(n), sorted_data(n))
        allocate (stack_left(64), stack_right(64))
        allocate (dummy_perm(n))
        ! For reproducibility: initialize the random number generator seed
        call random_seed(size=n_seed)
        allocate (seed_array(n_seed))
        seed_array = 42  ! Fixed value for reproducibility
        call random_seed(put=seed_array)
        deallocate (seed_array)
        call random_number(rdata)
        data = int(rdata*10000)
        ! Create reference sorted array using built-in approach
        sorted = data
        dummy_perm = [(i, i=1, n)]
        call sort_array(sorted, dummy_perm, stack_left, stack_right)
        sorted = sorted(dummy_perm)
        ! Test our sorting
        perm = [(i, i=1, n)]
        call sort_array(data, perm, stack_left, stack_right)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, sorted, n, "test_sort_large_random: sorted values mismatch")
        deallocate (rdata, data, perm, sorted, sorted_data, stack_left, stack_right, dummy_perm)
    end subroutine test_sort_large_random

    !> Test heapsort implementation for real arrays (uses f42_sort_impl heapsort)
    subroutine test_heapsort_real()
        real(real64), dimension(5) :: data = [5.0d0, 3.0d0, 4.0d0, 1.0d0, 2.0d0]
        integer(int32), dimension(5) :: perm
        real(real64), dimension(5) :: expected = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
        real(real64), dimension(5) :: sorted_data
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_real_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected, 5, 1d-12, "test_heapsort_real: values not sorted")
    end subroutine test_heapsort_real

    !> Test heapsort implementation for integer arrays
    subroutine test_heapsort_integer()
        integer(int32), dimension(5) :: data = [5, 3, 4, 1, 2]
        integer(int32), dimension(5) :: perm, expected = [1, 2, 3, 4, 5]
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected, 5, "test_heapsort_integer: values not sorted")
    end subroutine test_heapsort_integer

    !> Test heapsort implementation for character arrays
    subroutine test_heapsort_character()
        character(len=1), dimension(5) :: data = ["d", "b", "c", "a", "e"]
        integer(int32), dimension(5) :: perm
        character(len=1), dimension(5) :: expected = ["a", "b", "c", "d", "e"]
        character(len=1), dimension(5) :: sorted_data
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_character_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_char(sorted_data, expected, 1, 5, "test_heapsort_character: values not sorted")
    end subroutine test_heapsort_character

    !> Heapsort variant: integer ascending (mirrors test_sort_integer_ascending)
    subroutine test_heapsort_integer_ascending()
        integer(int32), dimension(5) :: data = [5, 2, 9, 1, 6]
        integer(int32), dimension(5) :: perm, expected_sorted = [1, 2, 5, 6, 9]
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 5, "test_heapsort_integer_ascending: sorted values mismatch")
    end subroutine test_heapsort_integer_ascending

    !> Heapsort variant: real descending (mirrors test_sort_real_descending)
    subroutine test_heapsort_real_descending()
        real(real64), dimension(5) :: data = [3.5d0, 2.2d0, 8.8d0, 1.1d0, 7.7d0]
        real(real64), dimension(5) :: expected_sorted = [1.1d0, 2.2d0, 3.5d0, 7.7d0, 8.8d0]
        real(real64), dimension(5) :: sorted_data
        integer(int32), dimension(5) :: perm
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_real_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 5, 1d-12, "test_heapsort_real_descending: sorted values mismatch")
    end subroutine test_heapsort_real_descending

    !> Heapsort variant: character random (mirrors test_sort_char_random)
    subroutine test_heapsort_char_random()
        character(len=8), dimension(5) :: data = ['dog     ', 'apple   ', 'zebra   ', 'cat     ', 'bird    ']
        character(len=8), dimension(5) :: expected = ['apple   ', 'bird    ', 'cat     ', 'dog     ', 'zebra   ']
        integer(int32), dimension(5) :: perm
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_character_heapsort(data, perm)
        call assert_true(all(data(perm) == expected), "test_heapsort_char_random: sorted values mismatch")
    end subroutine test_heapsort_char_random

    !> Heapsort variant: sorted stability (mirrors test_sort_sorted_stability)
    subroutine test_heapsort_sorted_stability()
        integer(int32), dimension(5) :: data = [1, 2, 3, 4, 5]
        integer(int32), dimension(5) :: expected = [1, 2, 3, 4, 5]
        integer(int32), dimension(5) :: perm
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected, 5, "test_heapsort_sorted_stability: sorted values mismatch")
    end subroutine test_heapsort_sorted_stability

    !> Heapsort variant: empty array (mirrors test_sort_empty_array)
    subroutine test_heapsort_empty_array()
        integer(int32), dimension(0) :: data
        integer(int32), dimension(0) :: perm

        ! Should not crash
        call sort_integer_heapsort(data, perm)
        call assert_true(.true., "test_heapsort_empty_array: empty array handling")
    end subroutine test_heapsort_empty_array

    !> Heapsort variant: large random (mirrors test_sort_large_random)
    subroutine test_heapsort_large_random()
        real(real64), allocatable :: rdata(:)
        integer(int32), allocatable :: data(:), perm(:), sorted(:), sorted_data(:)
        integer(int32) :: n
        integer(int32), allocatable :: stack_left(:), stack_right(:), dummy_perm(:)
        integer(int32) :: i
        integer(int32) :: n_seed
        integer(int32), allocatable :: seed_array(:)
        n = 1000
        allocate (rdata(n), data(n), perm(n), sorted(n), sorted_data(n))
        allocate (dummy_perm(n))
        allocate (stack_left(64), stack_right(64))
        ! For reproducibility: initialize the random number generator seed
        call random_seed(size=n_seed)
        allocate (seed_array(n_seed))
        seed_array = 42
        call random_seed(put=seed_array)
        deallocate (seed_array)
        call random_number(rdata)
        data = int(rdata*10000)
        sorted = data
        dummy_perm = [(i, i=1, n)]
        call sort_array(sorted, dummy_perm, stack_left, stack_right)
        sorted = sorted(dummy_perm)
        perm = [(i, i=1, n)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, sorted, n, "test_heapsort_large_random: sorted values mismatch")
        deallocate (rdata, data, perm, sorted, sorted_data, stack_left, stack_right, dummy_perm)
    end subroutine test_heapsort_large_random

    !> Test sorting when the input contains duplicate values for integers.
    subroutine test_sort_duplicates()
        integer(int32), dimension(6) :: data = [2, 1, 2, 4, 2, 3]
        integer(int32), dimension(6) :: perm, expected_sorted = [1, 2, 2, 2, 3, 4]
        integer(int32), dimension(6) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 6)]
        call sort_array(data, perm, stack_left, stack_right)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 6, "test_sort_duplicates (quicksort): sorted values mismatch")

        ! Also exercise heapsort implementation on the same data
        perm = [(i, i=1, 6)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 6, "test_sort_duplicates (heapsort): sorted values mismatch")
    end subroutine test_sort_duplicates

    !> Test sorting with negative numbers for integers.
    subroutine test_sort_negatives()
        integer(int32), dimension(5) :: data = [3, -1, 5, 0, -10]
        integer(int32), dimension(5) :: perm, expected_sorted = [-10, -1, 0, 3, 5]
        integer(int32), dimension(5) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 5, "test_sort_negatives (quicksort): sorted values mismatch")

        perm = [(i, i=1, 5)]
        call sort_integer_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_int(sorted_data, expected_sorted, 5, "test_sort_negatives (heapsort): sorted values mismatch")
    end subroutine test_sort_negatives

    !> Test sorting when the input contains duplicate values for reals.
    subroutine test_sort_duplicates_real()
        real(real64), dimension(6) :: data = [2.0d0, 1.0d0, 2.0d0, 4.0d0, 2.0d0, 3.0d0]
        integer(int32), dimension(6) :: perm
        real(real64), dimension(6) :: expected_sorted = [1.0d0, 2.0d0, 2.0d0, 2.0d0, 3.0d0, 4.0d0]
        real(real64), dimension(6) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 6)]
        call sort_array(data, perm, stack_left, stack_right)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 6, 1d-12, &
                                     "test_sort_duplicates_real (quicksort): sorted values mismatch")

        perm = [(i, i=1, 6)]
        call sort_real_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 6, 1d-12, &
                                     "test_sort_duplicates_real (heapsort): sorted values mismatch")
    end subroutine test_sort_duplicates_real

    !> Test sorting with negative numbers for reals.
    subroutine test_sort_negatives_real()
        real(real64), dimension(5) :: data = [3.0d0, -1.0d0, 5.0d0, 0.0d0, -10.0d0]
        integer(int32), dimension(5) :: perm
        real(real64), dimension(5) :: expected_sorted = [-10.0d0, -1.0d0, 0.0d0, 3.0d0, 5.0d0]
        real(real64), dimension(5) :: sorted_data
        integer(int32) :: stack_left(20), stack_right(20)
        integer(int32) :: i

        perm = [(i, i=1, 5)]
        call sort_array(data, perm, stack_left, stack_right)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 5, 1d-12, &
                                     "test_sort_negatives_real (quicksort): sorted values mismatch")

        perm = [(i, i=1, 5)]
        call sort_real_heapsort(data, perm)
        sorted_data = data(perm)
        call assert_equal_array_real(sorted_data, expected_sorted, 5, 1d-12, &
                                     "test_sort_negatives_real (heapsort): sorted values mismatch")
    end subroutine test_sort_negatives_real

    !> Test behavior when the array contains a NaN value: NaN should end up as the last element.
    subroutine test_sort_nan()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_is_nan, ieee_quiet_nan
        real(real64), dimension(4) :: data
        real(real64), dimension(3) :: expected_non_nan = [1.0d0, 2.0d0, 3.0d0]
        real(real64), dimension(3) :: sorted_non_nan
        integer(int32), dimension(4) :: perm
        integer(int32) :: stack_left(20), stack_right(20)
        real(real64) :: nanval
        integer(int32) :: n
        integer(int32) :: i

        nanval = ieee_value(0.0_real64, ieee_quiet_nan)
        data = [3.0d0, 1.0d0, nanval, 2.0d0]
        n = size(data)

        ! Quicksort variant
        perm = [(i, i=1, n)]
        call sort_array(data, perm, stack_left, stack_right)
        call assert_true(ieee_is_nan(data(perm(n))), "test_sort_nan (quicksort): NaN not last")
        sorted_non_nan = data(perm(1:n - 1))
        call assert_equal_array_real(sorted_non_nan, expected_non_nan, n - 1, 1d-12, &
                                     "test_sort_nan (quicksort): non-NaN values not sorted")

        ! Heapsort variant
        data = [3.0d0, 1.0d0, nanval, 2.0d0]
        perm = [(i, i=1, n)]
        call sort_real_heapsort(data, perm)
        call assert_true(ieee_is_nan(data(perm(n))), "test_sort_nan (heapsort): NaN not last")
        sorted_non_nan = data(perm(1:n - 1))
        call assert_equal_array_real(sorted_non_nan, expected_non_nan, n - 1, 1d-12, &
                                     "test_sort_nan (heapsort): non-NaN values not sorted")
    end subroutine test_sort_nan

    !> top_k_heap_push on a 3-slot heap seeded with 0, values 5, 1, -1, 7, 3, 3 pushed in turn,
    !| every heap state derived by hand:
    !| - 5 displaces the root 0: both children 0, the left is taken -> [0, 5, 0];
    !| - 1 displaces the root 0: the right child 0 is smaller than the left 5 -> [0, 5, 1];
    !| - -1 is not above the root 0: dropped, heap unchanged;
    !| - 7 displaces the root 0: the right child 1 is the smaller -> [1, 5, 7];
    !| - 3 displaces the root 1: the smaller child 5 is not below 3, the sift stops at once
    !|   -> [3, 5, 7];
    !| - 3 again equals the root: not above it, dropped.
    !| The heap then holds the 3 largest values pushed {3, 5, 7}, its root the smallest of them.
    !| A 2-slot heap exercises the left-child-only sift: [0, 0] + 4 -> [0, 4], + 2 -> [2, 4].
    subroutine test_top_k_heap_push()
        character(len=*), parameter :: name = "test_top_k_heap_push: "
        real(real64) :: heap(3), heap2(2)

        heap = 0.0_real64
        call top_k_heap_push(heap, 3_int32, 5.0_real64)
        call assert_equal_array_real(heap, [0.0_real64, 5.0_real64, 0.0_real64], 3_int32, 0.0_real64, name//"push 5")
        call top_k_heap_push(heap, 3_int32, 1.0_real64)
        call assert_equal_array_real(heap, [0.0_real64, 5.0_real64, 1.0_real64], 3_int32, 0.0_real64, name//"push 1")
        call top_k_heap_push(heap, 3_int32, -1.0_real64)
        call assert_equal_array_real(heap, [0.0_real64, 5.0_real64, 1.0_real64], 3_int32, 0.0_real64, &
                                     name//"push -1 dropped")
        call top_k_heap_push(heap, 3_int32, 7.0_real64)
        call assert_equal_array_real(heap, [1.0_real64, 5.0_real64, 7.0_real64], 3_int32, 0.0_real64, name//"push 7")
        call top_k_heap_push(heap, 3_int32, 3.0_real64)
        call assert_equal_array_real(heap, [3.0_real64, 5.0_real64, 7.0_real64], 3_int32, 0.0_real64, &
                                     name//"push 3, sift stops at once")
        call top_k_heap_push(heap, 3_int32, 3.0_real64)
        call assert_equal_array_real(heap, [3.0_real64, 5.0_real64, 7.0_real64], 3_int32, 0.0_real64, &
                                     name//"push 3 equal to the root dropped")

        heap2 = 0.0_real64
        call top_k_heap_push(heap2, 2_int32, 4.0_real64)
        call assert_equal_array_real(heap2, [0.0_real64, 4.0_real64], 2_int32, 0.0_real64, name//"2 slots, push 4")
        call top_k_heap_push(heap2, 2_int32, 2.0_real64)
        call assert_equal_array_real(heap2, [2.0_real64, 4.0_real64], 2_int32, 0.0_real64, name//"2 slots, push 2")
    end subroutine test_top_k_heap_push

    !> bottom_k_heap_push, the mirror of test_top_k_heap_push: a 3-slot heap seeded with 10,
    !| values 5, 9, 11, 3, 7, 7 pushed in turn:
    !| - 5 displaces the root 10: both children 10, the left is taken -> [10, 5, 10];
    !| - 9 displaces the root 10: the right child 10 is larger than the left 5 -> [10, 5, 9];
    !| - 11 is not below the root 10: dropped;
    !| - 3 displaces the root 10: the right child 9 is the larger -> [9, 5, 3];
    !| - 7 displaces the root 9: the larger child 5 is not above 7, the sift stops -> [7, 5, 3];
    !| - 7 again equals the root: dropped.
    !| The heap then holds the 3 smallest values pushed {3, 5, 7}. A 2-slot heap exercises the
    !| left-child-only sift: [10, 10] + 4 -> [10, 4], + 6 -> [6, 4].
    subroutine test_bottom_k_heap_push()
        character(len=*), parameter :: name = "test_bottom_k_heap_push: "
        real(real64) :: heap(3), heap2(2)

        heap = 10.0_real64
        call bottom_k_heap_push(heap, 3_int32, 5.0_real64)
        call assert_equal_array_real(heap, [10.0_real64, 5.0_real64, 10.0_real64], 3_int32, 0.0_real64, name//"push 5")
        call bottom_k_heap_push(heap, 3_int32, 9.0_real64)
        call assert_equal_array_real(heap, [10.0_real64, 5.0_real64, 9.0_real64], 3_int32, 0.0_real64, name//"push 9")
        call bottom_k_heap_push(heap, 3_int32, 11.0_real64)
        call assert_equal_array_real(heap, [10.0_real64, 5.0_real64, 9.0_real64], 3_int32, 0.0_real64, &
                                     name//"push 11 dropped")
        call bottom_k_heap_push(heap, 3_int32, 3.0_real64)
        call assert_equal_array_real(heap, [9.0_real64, 5.0_real64, 3.0_real64], 3_int32, 0.0_real64, name//"push 3")
        call bottom_k_heap_push(heap, 3_int32, 7.0_real64)
        call assert_equal_array_real(heap, [7.0_real64, 5.0_real64, 3.0_real64], 3_int32, 0.0_real64, &
                                     name//"push 7, sift stops at once")
        call bottom_k_heap_push(heap, 3_int32, 7.0_real64)
        call assert_equal_array_real(heap, [7.0_real64, 5.0_real64, 3.0_real64], 3_int32, 0.0_real64, &
                                     name//"push 7 equal to the root dropped")

        heap2 = 10.0_real64
        call bottom_k_heap_push(heap2, 2_int32, 4.0_real64)
        call assert_equal_array_real(heap2, [10.0_real64, 4.0_real64], 2_int32, 0.0_real64, name//"2 slots, push 4")
        call bottom_k_heap_push(heap2, 2_int32, 6.0_real64)
        call assert_equal_array_real(heap2, [6.0_real64, 4.0_real64], 2_int32, 0.0_real64, name//"2 slots, push 6")
    end subroutine test_bottom_k_heap_push

    !> A heap size of 0 leaves the array untouched for both pushes (the early return); a single
    !| slot is replaced by a value beyond it and kept otherwise, the sift exiting at once with no
    !| child.
    subroutine test_k_heap_push_size_zero_and_one()
        character(len=*), parameter :: name = "test_k_heap_push_size_zero_and_one: "
        real(real64) :: slot(1)

        slot = 1.0_real64
        call top_k_heap_push(slot, 0_int32, 5.0_real64)
        call assert_equal_real(slot(1), 1.0_real64, 0.0_real64, name//"top-k, size 0 writes nothing")
        call bottom_k_heap_push(slot, 0_int32, -5.0_real64)
        call assert_equal_real(slot(1), 1.0_real64, 0.0_real64, name//"bottom-k, size 0 writes nothing")

        call top_k_heap_push(slot, 1_int32, 5.0_real64)
        call assert_equal_real(slot(1), 5.0_real64, 0.0_real64, name//"top-k, size 1 takes a larger value")
        call top_k_heap_push(slot, 1_int32, 2.0_real64)
        call assert_equal_real(slot(1), 5.0_real64, 0.0_real64, name//"top-k, size 1 keeps its value over a smaller")
        call bottom_k_heap_push(slot, 1_int32, -5.0_real64)
        call assert_equal_real(slot(1), -5.0_real64, 0.0_real64, name//"bottom-k, size 1 takes a smaller value")
        call bottom_k_heap_push(slot, 1_int32, 0.0_real64)
        call assert_equal_real(slot(1), -5.0_real64, 0.0_real64, name//"bottom-k, size 1 keeps its value over a larger")
    end subroutine test_k_heap_push_size_zero_and_one

    !> The two sifts on a 5-element heap whose root alone breaks the heap condition: the root
    !| travels down two levels, each time to the child that keeps the heap valid.
    !| Min-heap [9, 1, 2, 3, 4]: 9 swaps with 1 (left), then with 3 (left of 2) -> [1, 3, 2, 9, 4].
    !| Max-heap [0, 8, 9, 7, 6]: 0 swaps with 9 (right), then stops (position 3 has no child)
    !| -> [9, 8, 0, 7, 6].
    subroutine test_heap_sift_down()
        character(len=*), parameter :: name = "test_heap_sift_down: "
        real(real64) :: heap(5)

        heap = [9.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        call minheap_sift_down(heap, 5_int32)
        call assert_equal_array_real(heap, [1.0_real64, 3.0_real64, 2.0_real64, 9.0_real64, 4.0_real64], 5_int32, &
                                     0.0_real64, name//"min-heap")

        heap = [0.0_real64, 8.0_real64, 9.0_real64, 7.0_real64, 6.0_real64]
        call maxheap_sift_down(heap, 5_int32)
        call assert_equal_array_real(heap, [9.0_real64, 8.0_real64, 0.0_real64, 7.0_real64, 6.0_real64], 5_int32, &
                                     0.0_real64, name//"max-heap")
    end subroutine test_heap_sift_down

end module mod_test_sorting
