! filepath: test/mod_test_kd_tree.f90
!> Unit test suite for kd_tree module.
module mod_test_kd_tree
  use kd_tree
  use f42_utils
  use tox_errors, only: set_ok, is_ok
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
    type(test_case) :: all_tests(12)
    
    all_tests(1) = test_case("test_kd_2d_cartesian", test_kd_2d_cartesian)
    all_tests(2) = test_case("test_kd_3d_spherical", test_kd_3d_spherical)
    all_tests(3) = test_case("test_kd_empty_array", test_kd_empty_array)
    all_tests(4) = test_case("test_kd_single_point", test_kd_single_point)
    all_tests(5) = test_case("test_kd_identical_points", test_kd_identical_points)
    all_tests(6) = test_case("test_kd_unit_vectors", test_kd_unit_vectors)
    all_tests(7) = test_case("test_kd_high_dim_low_points", test_kd_high_dim_low_points)
    all_tests(8) = test_case("test_kd_1d_sorted", test_kd_1d_sorted)
    all_tests(9) = test_case("test_kd_2d_minimal", test_kd_2d_minimal)
    all_tests(10) = test_case("test_kd_1d_minimal", test_kd_1d_minimal)
    all_tests(11) = test_case("test_kd_3d_large", test_kd_3d_large)
    all_tests(12) = test_case("test_kd_5d_medium", test_kd_5d_medium)
  end function get_all_tests

  !> Run all KD-Tree tests.
  subroutine run_all_tests_kd_tree()
    type(test_case) :: all_tests(12)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All KD-Tree tests passed successfully."
  end subroutine run_all_tests_kd_tree

  !> Run specific KD-Tree tests by name.
  subroutine run_named_tests_kd_tree(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(12)
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
  end subroutine run_named_tests_kd_tree

  !> Test 2D Cartesian KD-Tree.
  subroutine test_kd_2d_cartesian()
    integer(int32), parameter :: d = 2, n = 6
    real(real64) :: X(d,n) = reshape([1.0d0, 2.0d0, 2.0d0, 3.0d0, 3.0d0, 1.0d0, &
                                    4.0d0, 0.0d0, 0.0d0, 4.0d0, 5.0d0, 2.0d0], [d, n])
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for 2D cartesian tree: ', ierr
      error stop
    end if
    
    call assert_permutation(kd_ix, n, "2D Cartesian KD-Tree")
  end subroutine test_kd_2d_cartesian

  !> Test 3D Spherical KD-Tree.
  subroutine test_kd_3d_spherical()
    integer(int32), parameter :: d = 3, n = 8
    real(real64) :: V(d,n)
    integer(int32) :: sphere_ix(n), dim_order(d) = [1, 2, 3]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    call random_unit_vectors(V, d, n)
    call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build spherical kd failed for 3D spherical tree: ', ierr
      error stop
    end if
    call assert_permutation(sphere_ix, n, "3D Spherical KD-Tree")
  end subroutine test_kd_3d_spherical

  !> Test KD-Tree with empty array.
  subroutine test_kd_empty_array()
    integer(int32), parameter :: d = 2, n = 0
    real(real64) :: X(d,n)
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    ! This should return ERR_EMPTY_INPUT
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for empty input - should return ERR_EMPTY_INPUT but got ERR_OK'
      error stop
    end if
    
    call assert_true(.true., "KD-Tree empty array handling")
  end subroutine test_kd_empty_array

  !> Test KD-Tree with single point.
  subroutine test_kd_single_point()
    integer(int32), parameter :: d = 2, n = 1
    real(real64) :: X(d,n) = reshape([1.0d0, 2.0d0], [d, n])
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for single point input: ', ierr
      error stop
    end if
    
    call assert_equal_int(kd_ix(1), 1, "KD-Tree single point index incorrect")
  end subroutine test_kd_single_point

  !> Test KD-Tree with identical points.
  subroutine test_kd_identical_points()
    integer(int32), parameter :: d = 3, n = 5
    real(real64) :: X(d,n) = 1.0d0
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for identical points: ', ierr
      error stop
    end if

    call assert_permutation(kd_ix, n, "KD-Tree identical points")
  end subroutine test_kd_identical_points

  !> Test KD-Tree with unit vectors.
  subroutine test_kd_unit_vectors()
    integer(int32), parameter :: d = 4, n = 4
    real(real64) :: V(d,n) = 0.0d0
    integer(int32) :: sphere_ix(n), dim_order(d) = [1, 2, 3, 4]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n)
    real(real64) :: subarray(n)
    integer(int32) :: i, ierr
    integer(int32) :: recursion_stack(3, n)

    call set_ok(ierr)
    
    do i = 1, d
      V(i,i) = 1.0d0
    end do
    
    call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build shperical kd failed for unit vectors: ', ierr
      error stop ierr
    end if

    call assert_permutation(sphere_ix, n, "KD-Tree unit vectors")
  end subroutine test_kd_unit_vectors

  !> Test KD-Tree with high dimension and few points.
  subroutine test_kd_high_dim_low_points()
    integer(int32), parameter :: d = 10, n = 3
    real(real64) :: X(d,n)
    integer :: i
    integer(int32) :: kd_ix(n), dim_order(d) = [(i, i = 1, d)]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3, n)
    call set_ok(ierr)
    
    call random_matrix(X, d, n)
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for high dimensioninality: ', ierr
      error stop 
    end if
    call assert_permutation(kd_ix, n, "KD-Tree high dimension low points")
  end subroutine test_kd_high_dim_low_points

  !> Test 1D sorted KD-Tree.
  subroutine test_kd_1d_sorted()
    integer(int32), parameter :: d = 1, n = 10
    real(real64) :: X(d,n)
    integer(int32) :: kd_ix(n), dim_order(d) = [1]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n)
    real(real64) :: subarray(n)
    integer(int32) :: i, ierr
    integer(int32) :: recursion_stack(3, n)

    call set_ok(ierr)
    
    do i = 1, n
      X(1,i) = i
    end do
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for sorted 1D array: ', ierr
      error stop
    end if

    call assert_permutation(kd_ix, n, "1D sorted KD-Tree")
  end subroutine test_kd_1d_sorted

  !> Test 2D minimal KD-Tree.
  subroutine test_kd_2d_minimal()
    integer(int32), parameter :: d = 2, n = 2
    real(real64) :: X(d,n) = reshape([1.0d0, 2.0d0, 2.0d0, 1.0d0], [d, n])
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3,n)

    call set_ok(ierr)
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for 2D minimal test: ', ierr
      error stop 
    end if

    call assert_permutation(kd_ix, n, "2D minimal KD-Tree")
  end subroutine test_kd_2d_minimal

  !> Test 1D minimal KD-Tree.
  subroutine test_kd_1d_minimal()
    integer(int32), parameter :: d = 1, n = 2
    real(real64) :: X(d,n) = reshape([1.0d0, 2.0d0], [d, n])
    integer(int32) :: kd_ix(n), dim_order(d) = [1]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3, n)

    call set_ok(ierr)
    
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for 1D minimal test: ', ierr
      error stop
    end if

    call assert_permutation(kd_ix, n, "1D minimal KD-Tree")
  end subroutine test_kd_1d_minimal

  !> Test 3D large KD-Tree.
  subroutine test_kd_3d_large()
    integer(int32), parameter :: d = 3, n = 100
    real(real64) :: X(d,n)
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    real(real64) :: val(d)
    integer(int32) :: recursion_stack(3, n)

    call set_ok(ierr)
    
    call random_matrix(X, d, n)
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for 3D large input: ', ierr
      error stop
    end if
    call get_kd_point(X, kd_ix, 4, val, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Get kd point failed for 3D large input: ', ierr
      error stop
    end if

    call assert_permutation(kd_ix, n, "3D large KD-Tree")
  end subroutine test_kd_3d_large

  !> Test 5D medium KD-Tree.
  subroutine test_kd_5d_medium()
    integer(int32), parameter :: d = 5, n = 10
    real(real64) :: X(d,n)
    integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3, 4, 5]
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), ierr
    real(real64) :: subarray(n)
    integer(int32) :: recursion_stack(3, n)

    call set_ok(ierr)
    
    call random_matrix(X, d, n)
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right, recursion_stack, ierr)
    if(.not. is_ok(ierr)) then 
      write(*,*) 'Build kd index failed for 5D tree: ', ierr
      error stop
    end if
    call assert_permutation(kd_ix, n, "5D medium KD-Tree")
  end subroutine test_kd_5d_medium

  !> Helper: Generate random unit vectors.
  subroutine random_unit_vectors(V, d, n)
    integer(int32), intent(in) :: d, n
    real(real64), intent(out) :: V(d,n)
    integer(int32) :: i
    real(real64) :: norm
    call random_seed()
    do i = 1, n
      call random_number(V(:,i))
      V(:,i) = V(:,i) - 0.5d0
      norm = sqrt(sum(V(:,i)**2))
      if (norm > 0) V(:,i) = V(:,i) / norm
    end do
  end subroutine random_unit_vectors

  !> Helper: Generate random matrix.
  subroutine random_matrix(X, d, n)
    integer(int32), intent(in) :: d, n
    real(real64), intent(out) :: X(d,n)
    integer(int32) :: j
    call random_seed()
    do j = 1, n
      call random_number(X(:,j))
    end do
  end subroutine random_matrix

  !> Helper: Convert integer to string.
  function str(i) result(s)
    integer(int32), intent(in) :: i
    character(len=32) :: s
    write(s, *) i
    s = adjustl(s)
  end function str

end module mod_test_kd_tree