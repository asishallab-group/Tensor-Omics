! filepath: test/mod_test_bst.f90
!> Unit test suite for binary_search_tree module.
module mod_test_bst
  use binary_search_tree
  use tox_errors, only: is_ok, set_ok
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
    type(test_case) :: all_tests(7)
    
    all_tests(1) = test_case("test_bst_index_construction", test_bst_index_construction)
    all_tests(2) = test_case("test_bst_sorted_values", test_bst_sorted_values)
    all_tests(3) = test_case("test_bst_range_query", test_bst_range_query)
    all_tests(4) = test_case("test_bst_empty_array", test_bst_empty_array)
    all_tests(5) = test_case("test_bst_single_element", test_bst_single_element)
    all_tests(6) = test_case("test_bst_identical_values", test_bst_identical_values)
    all_tests(7) = test_case("test_bst_large_random", test_bst_large_random)
  end function get_all_tests

  !> Run all BST tests.
  subroutine run_all_tests_bst()
    type(test_case) :: all_tests(7)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All BST tests passed successfully."
  end subroutine run_all_tests_bst

  !> Run specific BST tests by name.
  subroutine run_named_tests_bst(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(7)
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
  end subroutine run_named_tests_bst

  !> Test BST index construction and monotonicity.
  subroutine test_bst_index_construction()
    integer(int32), parameter :: n = 100
    real(real64) :: x(n)
    integer(int32) :: ix(n), stack_left(n), stack_right(n)
    integer(int32) :: i, ierr
    logical :: is_sorted

    call set_ok(ierr)
    
    call random_array(x, n)
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) error stop
    ! Check monotonicity of x(ix)
    is_sorted = .true.
    do i = 2, n
      if (x(ix(i)) < x(ix(i-1))) then
        is_sorted = .false.
        exit
      end if
    end do
    
    call assert_true(is_sorted, "BST index test FAILED: x(ix) is not monotonic.")
  end subroutine test_bst_index_construction

  !> Test get_sorted_value function.
  subroutine test_bst_sorted_values()
    integer(int32), parameter :: n = 10
    real(real64) :: x(n) = [3.0d0, 1.0d0, 4.0d0, 2.0d0, 5.0d0, 7.0d0, 6.0d0, 9.0d0, 8.0d0, 10.0d0]
    integer(int32) :: ix(n), stack_left(n), stack_right(n)
    real(real64) :: val
    integer(int32) :: ierr

    call set_ok(ierr)
    
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Build bst index failed for sorted values: ', ierr
      error stop 
    end if 
    val = get_sorted_value(x, ix, 3, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Get sorted value failed: ', ierr
      error stop
    end if
    
    call assert_equal_real(val, 3.0d0, 1d-12, "get_sorted_value returned incorrect value")
  end subroutine test_bst_sorted_values

  !> Test BST range query functionality.
  subroutine test_bst_range_query()
    integer(int32), parameter :: n = 10
    real(real64) :: x(n) = [3.0d0, 1.0d0, 4.0d0, 2.0d0, 5.0d0, 7.0d0, 6.0d0, 9.0d0, 8.0d0, 10.0d0]
    integer(int32) :: ix(n), stack_left(n), stack_right(n)
    integer(int32) :: res_ix(n), res_n, ierr

    call set_ok(ierr)
    
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build bst index failed: ', ierr
      error stop
    end if
    call bst_range_query(x, ix, n, 2.5d0, 7.5d0, res_ix, res_n, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Bst range query failed: ', ierr
      error stop
    end if

    call assert_true(res_n == 5, "BST range query returned incorrect count")
  end subroutine test_bst_range_query

  !> Test BST with empty array.
  subroutine test_bst_empty_array()
    integer(int32), parameter :: n = 0
    real(real64) :: x(n)
    integer(int32) :: ix(n), stack_left(n), stack_right(n), ierr

    call set_ok(ierr)

    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(is_ok(ierr)) then 
      write(*,*) 'Build bst index failed for empty array: Expected ERR_EMPTY_INPUT but got ERR_OK '
      error stop
    end if
    call assert_true(.true., "BST empty array handling")
  end subroutine test_bst_empty_array

  !> Test BST with single element.
  subroutine test_bst_single_element()
    integer(int32), parameter :: n = 1
    real(real64) :: x(n) = [42.0d0]
    integer(int32) :: ix(n), stack_left(n), stack_right(n), ierr
    
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Build bst index failed for single element: ', ierr
      error stop
    end if
    call assert_equal_int(ix(1), 1, "BST single element index incorrect")
    call assert_equal_real(x(ix(1)), 42.0d0, 1d-12, "BST single element value incorrect")
  end subroutine test_bst_single_element

  !> Test BST with identical values.
  subroutine test_bst_identical_values()
    integer(int32), parameter :: n = 5
    real(real64) :: x(n) = 7.0d0
    integer(int32) :: ix(n), stack_left(n), stack_right(n)
    integer(int32) :: i, ierr

    call set_ok(ierr)
    
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build bst index failed for identical values: ', ierr
      error stop
    end if
    ! Should still be a valid permutation
    do i = 1, n
      call assert_true(ix(i) >= 1 .and. ix(i) <= n, "BST identical values index out of bounds")
    end do
  end subroutine test_bst_identical_values

  !> Test BST with large random array.
  subroutine test_bst_large_random()
    integer(int32), parameter :: n = 1000
    real(real64) :: x(n)
    integer(int32) :: ix(n), stack_left(n), stack_right(n)
    integer(int32) :: i, ierr
    logical :: is_sorted
    
    call set_ok(ierr)
    call random_array(x, n)
    call build_bst_index(x, n, ix, stack_left, stack_right, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build bst index failed for large random values: ', ierr
      error stop    ! Should still be a valid permutation
    end if
    ! Check monotonicity of x(ix)
    is_sorted = .true.
    do i = 2, n
      if (x(ix(i)) < x(ix(i-1))) then
        is_sorted = .false.
        exit
      end if
    end do
    
    call assert_true(is_sorted, "Large BST index test FAILED: x(ix) is not monotonic.")
  end subroutine test_bst_large_random

  !> Fill an array with random real values in [0,1).
  subroutine random_array(arr, nval)
    real(real64), intent(out) :: arr(:)
    integer(int32), intent(in) :: nval
    integer(int32) :: j
    call random_seed()
    do j = 1, nval
      call random_number(arr(j))
    end do
  end subroutine random_array

end module mod_test_bst