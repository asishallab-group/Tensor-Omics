! filepath: test/mod_test_sorting.f90
!> Unit test suite for f42_utils module.
module mod_test_sorting
  use f42_utils
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(9)
    
    all_tests(1) = test_case("test_sort_real", test_sort_real)
    all_tests(2) = test_case("test_sort_integer", test_sort_integer)
    all_tests(3) = test_case("test_sort_character", test_sort_character)
    all_tests(4) = test_case("test_sort_integer_ascending", test_sort_integer_ascending)
    all_tests(5) = test_case("test_sort_real_descending", test_sort_real_descending)
    all_tests(6) = test_case("test_sort_char_random", test_sort_char_random)
    all_tests(7) = test_case("test_sort_sorted_stability", test_sort_sorted_stability)
    all_tests(8) = test_case("test_sort_empty_array", test_sort_empty_array)
    all_tests(9) = test_case("test_sort_large_random", test_sort_large_random)
  end function get_all_tests

  !> Run all sorting tests.
  subroutine run_all_tests_sorting()
    type(test_case) :: all_tests(9)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All sorting tests passed successfully."
  end subroutine run_all_tests_sorting

  !> Run specific sorting tests by name.
  subroutine run_named_tests_sorting(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(9)
    integer(int32) :: i, j
    logical :: found
    
    all_tests = get_all_tests()
    
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_sorting

  !> Test sorting of a real array using permutation vector.
  subroutine test_sort_real()
    real(real64), dimension(5) :: data = [3.0d0, 1.0d0, 5.0d0, 2.0d0, 4.0d0]
    integer(int32), dimension(5) :: perm, expected = [2, 4, 1, 5, 3]
    integer(int32) :: stack_left(20), stack_right(20), i
    
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_int(perm, expected, 5, "test_sort_real: permutation mismatch")
  end subroutine test_sort_real

  !> Test sorting of an integer array using permutation vector.
  subroutine test_sort_integer()
    integer(int32), dimension(4) :: data = [10, 3, 7, 1]
    integer(int32), dimension(4) :: perm, expected = [4, 2, 3, 1]
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 4)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_int(perm, expected, 4, "test_sort_integer: permutation mismatch")
  end subroutine test_sort_integer

  !> Test sorting of a character array using permutation vector.
  subroutine test_sort_character()
    character(len=6), dimension(3) :: data = ['delta ', 'alpha ', 'beta  ']
    integer(int32), dimension(3) :: perm, expected = [2, 3, 1]
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 3)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_int(perm, expected, 3, "test_sort_character: permutation mismatch")
  end subroutine test_sort_character

  !> Test that sorted values are in ascending order for integers.
  subroutine test_sort_integer_ascending()
    integer(int32), dimension(5) :: data = [5, 2, 9, 1, 6]
    integer(int32), dimension(5) :: perm, expected_sorted = [1, 2, 5, 6, 9]
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_int(data(perm), expected_sorted, 5, "test_sort_integer_ascending: sorted values mismatch")
  end subroutine test_sort_integer_ascending

  !> Test that sorted values are in ascending order for reals.
  subroutine test_sort_real_descending()
    real(real64), dimension(5) :: data = [3.5d0, 2.2d0, 8.8d0, 1.1d0, 7.7d0]
    real(real64), dimension(5) :: expected_sorted = [1.1d0, 2.2d0, 3.5d0, 7.7d0, 8.8d0]
    integer(int32), dimension(5) :: perm
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_real(data(perm), expected_sorted, 5, 1d-12, &
                            "test_sort_real_descending: sorted values mismatch")
  end subroutine test_sort_real_descending

  !> Test sorting of character array with lexicographic ordering.
  subroutine test_sort_char_random()
    character(len=8), dimension(5) :: data = ['dog     ', 'apple   ', 'zebra   ', 'cat     ', 'bird    ']
    character(len=8), dimension(5) :: expected = ['apple   ', 'bird    ', 'cat     ', 'dog     ', 'zebra   ']
    integer(int32), dimension(5) :: perm
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_true(all(data(perm) == expected), "test_sort_char_random: sorted values mismatch")
  end subroutine test_sort_char_random

  !> Test sorting stability with already sorted array.
  subroutine test_sort_sorted_stability()
    integer(int32), dimension(5) :: data = [1, 2, 3, 4, 5]
    integer(int32), dimension(5) :: expected = [1, 2, 3, 4, 5]
    integer(int32), dimension(5) :: perm
    integer(int32) :: stack_left(20), stack_right(20)
    integer(int32) :: i
    
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    
    call assert_equal_array_int(data(perm), expected, 5, "test_sort_sorted_stability: sorted values mismatch")
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
    integer(int32), allocatable :: data(:), perm(:), sorted(:)
    integer(int32) :: n
    integer(int32), allocatable :: stack_left(:), stack_right(:), dummy_perm(:)
    integer(int32) :: i
    integer(int32) :: n_seed
    integer(int32), allocatable :: seed_array(:)
    n = 1000
    allocate(rdata(n), data(n), perm(n), sorted(n))
    allocate(stack_left(64), stack_right(64))
    allocate(dummy_perm(n))
    ! For reproducibility: initialize the random number generator seed
    call random_seed(size=n_seed)
    allocate(seed_array(n_seed))
    seed_array = 42  ! Fixed value for reproducibility
    call random_seed(put=seed_array)
    deallocate(seed_array)
    call random_number(rdata)
    data = int(rdata * 10000)
    ! Create reference sorted array using built-in approach
    sorted = data
    dummy_perm = [(i, i = 1, n)]
    call sort_array(sorted, dummy_perm, stack_left, stack_right)
    sorted = sorted(dummy_perm)
    ! Test our sorting
    perm = [(i, i = 1, n)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(data(perm), sorted, n, "test_sort_large_random: sorted values mismatch")
    deallocate(rdata, data, perm, sorted, stack_left, stack_right, dummy_perm)
  end subroutine test_sort_large_random

end module mod_test_sorting