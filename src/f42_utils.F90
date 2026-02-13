#include "macros.h"

!> Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.
module f42_utils
  use safeguard
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_DIVISION_BY_ZERO, set_ok, set_err_once, set_err, validate_in_range_real, is_err
  use, intrinsic :: ieee_arithmetic, only: ieee_next_after, ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_finite, ieee_is_nan
  implicit none

  public :: sort_real, sort_integer, sort_character
  public :: sort_array
  public :: compute_edf, compute_edf_alloc

  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array

  interface sort_array_heapsort
    module procedure sort_real_heapsort, sort_integer_heapsort, sort_character_heapsort
  end interface sort_array_heapsort

  real(real64), parameter :: PI = 4.0_real64 * atan(1.0_real64)
  real(real64), parameter :: EPS = epsilon(1.0_real64)
contains

  !> Compute logarithm for any base
  pure subroutine logx(val, base, exponent, ierr)
      real(real64), intent(in) :: val
        !! Value (`x` in $ b^y = x $)
      real(real64), intent(in) :: base
        !! Base (`b` in $ b^y = x $)
      real(real64), intent(out) :: exponent
        !! Exponent (`y` in $ b^y = x $)
      integer(int32), intent(out) :: ierr
        !! Error code

      call set_ok(ierr)
      call validate_in_range_real(val, ierr, min=above(0.0_real64))
      call validate_in_range_real(base, ierr, min=above(0.0_real64))
      if (is_close(base, 1.0_real64)) call set_err(ierr, ERR_DIVISION_BY_ZERO)
      if (is_err(ierr)) return
      
      exponent = log(val) / log(base)
  end subroutine logx

  !> Initialize Fortran's random number generator
  subroutine init_random(seed)
    integer(int32), intent(in), optional :: seed
      !! optional random seed, default: 42
      !!
      !! @note
      !! This subroutine uses the intrinsic `random_seed` subroutine that expects default kind integer.
      !! To map `int32` to default kind, it is being clamped by modulo to be in range.
      !! @endnote

    integer(int32) :: actual_seed

    ! IMPORTANT: these locals need to be default kind integer
    integer :: seed_default_kind, size, i
    integer(int32), parameter :: max_default_kind_val = int(huge(1), kind=int32)

    M_DEFAULT_VAL(seed, actual_seed, 42_int32)

    ! clamp to default kind range
    seed_default_kind = int(mod(actual_seed, max_default_kind_val))

    ! determine needed array size to seed random numbers
    call random_seed(size=size)

    ! create the seeding array, has negligible size
    call random_seed(put=[(seed_default_kind,i=1,size)])
  end subroutine init_random

  !> Returns a random real number `min <= rand_num < max`. If `min > max`, it will be `max <= rand_num < min`. If `min == max`, it will be `min`.
  real(real64) function rand_range(min, max) result(rand_num)
    real(real64), intent(in) :: min
    real(real64), intent(in) :: max

    call random_number(rand_num)
    rand_num = min + rand_num * (max - min)
  end function rand_range

  !> Returns the next representable float lower than a value. Helpful for exclusive upper bounds in ranges.
  pure real(real64) function below(val)
      real(real64), intent(in) :: val
  
      if (val == 0.0_real64) then
          below = -tiny(1.0_real64)
      else
          below = ieee_next_after(val, ieee_value(1.0_real64, ieee_negative_inf))
      end if
  end function below

  !> Returns the next representable float greater than a value. Helpful for exclusive upper bounds in ranges.
  pure real(real64) function above(val)
      real(real64), intent(in) :: val
  
      if (val == 0.0_real64) then
          above = tiny(1.0_real64)
      else
          above = ieee_next_after(val, ieee_value(1.0_real64, ieee_positive_inf))
      end if
  end function above

  pure logical function is_close(a, b)
    real(real64), intent(in) :: a
      !! First variable of comparison a==b
    real(real64), intent(in) :: b
      !! Second variable of comparison a==b

    real(real64) :: rel_tolerance

    if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
      rel_tolerance = EPS * max(abs(a), abs(b))
      is_close = abs(a - b) <= max(rel_tolerance, 1d-12)
    else
      is_close = a == b
    end if
  end function is_close

  !> Find the next power of two greater than or equal to n
  function next_power_of_two(n) result(power)
      integer(int32), intent(in) :: n
      !! input value
      integer(int32) :: power
      !! next greater value that is a power of two
      
      power = 2 ** (bit_size(n) - leadz(n - 1))
  end function next_power_of_two

  !> Computes the radian angle between two vectors
  pure subroutine angle_between(v1, v2, n_dims, angle, ierr)
    integer(int32), intent(in) :: n_dims
      !! number of elements in `v1` and `v2`
    real(real64), dimension(n_dims), intent(in) :: v1
      !! first vector for angle calculation
    real(real64), dimension(n_dims), intent(in) :: v2
      !! second vector for angle calculation
    real(real64), intent(out) :: angle
      !! will hold calculated angle
    integer(int32), intent(out) :: ierr
      !! Error code

    integer(int32) :: i
    real(real64) :: theta, dot_product, norm1, norm2, norm_product

    call set_ok(ierr)

    dot_product = 0
    norm1 = 0
    norm2 = 0
    do i = 1, size(v1)
      dot_product = dot_product + v1(i) * v2(i)
      norm1 = norm1 + v1(i) ** 2
      norm2 = norm2 + v2(i) ** 2
    end do
    norm_product = sqrt(norm1) * sqrt(norm2)
    if (is_close(norm_product, 0.0_real64)) then
        call set_err(ierr, ERR_DIVISION_BY_ZERO)
        return
    end if

    theta = dot_product / norm_product
    theta = max(-1.0_real64, min(1.0_real64, theta))
    angle = acos(theta)
  end subroutine angle_between

  !> Returns the given degrees in positive radian value (-90deg -> 3*PI/2, not -PI/2)
  pure real(real64) function radians(degrees)
    real(real64), intent(in) :: degrees
        !! degrees to be converted

    radians = modulo(degrees, 360.0_real64) * PI / 180 
  end function radians

  !> Calculates the euclidean norm of a vector
  pure real(real64) function norm(vector)
    real(real64), dimension(:), intent(in) :: vector
        !! Input vector the norm will be calcuated for

    integer(int32) :: i_dim

    norm = 0
    do i_dim = 1, size(vector)
      norm = norm + vector(i_dim) ** 2
    end do
    norm = sqrt(norm)
  end function norm

  !> Adds two vectors in-place
  pure subroutine add_vector(vector, to_be_added)
    real(real64), dimension(:), intent(inout) :: vector
        !! First vector, it will be modified in-place
    real(real64), dimension(:), intent(in) :: to_be_added
        !! Vector that should be added to `vector`

    integer(int32) :: i_dim

    do i_dim = 1, size(vector)
      vector(i_dim) = vector(i_dim) + to_be_added(i_dim)
    end do
  end subroutine add_vector

  !> Subtracts two vectors in-place
  pure subroutine subtract_vector(vector, to_be_subtracted)
    real(real64), dimension(:), intent(inout) :: vector
        !! First vector, it will be modified in-place
    real(real64), dimension(:), intent(in) :: to_be_subtracted
        !! Vector that should be subtracted from `vector`

    integer(int32) :: i_dim

    do i_dim = 1, size(vector)
      vector(i_dim) = vector(i_dim) - to_be_subtracted(i_dim)
    end do
  end subroutine subtract_vector
  
  !> Sort a real array indirectly using quicksort.
  !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
  pure subroutine sort_real(array, perm, stack_left, stack_right)
    !| Real input array to sort
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_real(array, perm, int(size(array), int32), stack_left, stack_right)
  end subroutine sort_real

  !> Sort an integer array indirectly using quicksort.
  !| Similar to `sort_real`, but for integer input.
  pure subroutine sort_integer(array, perm, stack_left, stack_right)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_int(array, perm, int(size(array), int32), stack_left, stack_right)
  end subroutine sort_integer

  !> Sort a character array indirectly using quicksort.
  !| Uses lexicographic ordering and permutation vector sorting.
  pure subroutine sort_character(array, perm, stack_left, stack_right)
    !| Character input array to sort
    character(len=*), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_char(array, perm, int(size(array), int32), stack_left, stack_right)
  end subroutine sort_character

  !> Sort a real array indirectly using heapsort.
  !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
  pure subroutine sort_real_heapsort(array, perm)
    !| Real input array to sort
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    call heapsort_real(array, perm)
  end subroutine sort_real_heapsort  

  !> Sort an integer array indirectly using heapsort.  
  !| Similar to `sort_real_heapsort`, but for integer input.
  pure subroutine sort_integer_heapsort(array, perm)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    call heapsort_integer(array, perm)
  end subroutine sort_integer_heapsort  

  !> Sort a character array indirectly using heapsort.  
  !| Uses lexicographic ordering and permutation vector sorting.  
  pure subroutine sort_character_heapsort(array, perm)
    !| Character input array to sort
    character(len=*), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    call heapsort_character(array, perm)
  end subroutine sort_character_heapsort  

  !> Internal quicksort implementation for real arrays.
  !| Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
  pure subroutine quicksort_real(array, perm, n, stack_left, stack_right)
    !| Real input array to sort
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)

    integer(int32) :: left, right, i, j, top, pivot_idx
    real(real64) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    ! Iterative quicksort using explicit stack
    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      ! Select pivot and initialize pointers
      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      ! Partitioning loop
      do
        do while (real_less(array(perm(i)), pivot_val))
          i = i + 1
        end do
        do while (real_greater(array(perm(j)), pivot_val))
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      ! Push new ranges onto stack
      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_real

  !> Internal quicksort implementation for integer arrays.
  !| Indirectly sorts `array` using `perm`, same algorithm as `quicksort_real`.
  pure subroutine quicksort_int(array, perm, n, stack_left, stack_right)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)

    integer(int32) :: left, right, i, j, top, pivot_idx
    integer(int32) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      do
        do while (array(perm(i)) < pivot_val)
          i = i + 1
        end do
        do while (array(perm(j)) > pivot_val)
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_int

  !> Internal quicksort implementation for character arrays.
  !| Lexicographic quicksort using string comparison, indirect via `perm`.
  pure subroutine quicksort_char(array, perm, n, stack_left, stack_right)
    !| Character input array to sort
    character(len=*), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    !| Temporary variables
    integer(int32) :: left, right, i, j, top, pivot_idx
    character(len=len(array)) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      do
        do while (array(perm(i)) < pivot_val)
          i = i + 1
        end do
        do while (array(perm(j)) > pivot_val)
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_char

  !Internal heapsort implementations for real arrays.
  !> Sorts indirectly using the permutation vector `perm`. Uses `heapify_real` to maintain heap property.
  pure subroutine heapsort_real(array, perm) 
    !| Real input array to sort     
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    integer(int32) :: n, i
    n = size(array)
    
    ! Build max heap
    do i = n / 2, 1, -1
      call heapify_real(array, perm, n, i)
    end do
    ! Heap sort
    do i = n, 2, -1
      call swap_int(perm(1), perm(i))
      call heapify_real(array, perm, i - 1, 1)
    end do

  contains

    !> Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_real(array, perm, heap_size, root)
      !| Real input array to sort
      real(real64), intent(in) :: array(:)
      !| Permutation vector that will be sorted
      integer(int32), intent(inout) :: perm(:)
      !| Size of the heap
      integer(int32), intent(in) :: heap_size
      !| Root index
      integer(int32), intent(in) :: root
      integer(int32) :: current, next_idx, largest_idx

      current = root

      do
        ! Use `largest_idx` directly as the left-child index: left = 2*current
        largest_idx = 2 * current
        ! If there is no left child the subtree is a leaf; we're done
        if (largest_idx > heap_size) exit

        next_idx = largest_idx + 1

        ! Only compare the right child (next_idx) when it actually exists
        if (next_idx <= heap_size) then
          if (real_greater(array(perm(next_idx)), array(perm(largest_idx)))) then
            largest_idx = next_idx
          end if
        end if

        ! If the larger child is greater than current, swap permutation entries
        if (real_greater(array(perm(largest_idx)), array(perm(current)))) then
          call swap_int(perm(current), perm(largest_idx))
          current = largest_idx
        else
          exit
        end if
      end do
    end subroutine heapify_real
  end subroutine heapsort_real

  !> Internal heapsort implementation for integer arrays.
  !| Indirectly sorts `array` using `perm`, same algorithm as `heapsort_real`.
  pure subroutine heapsort_integer(array, perm)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    integer(int32) :: n, i

    n = size(array)

    ! Build max-heap
    do i = n / 2, 1, -1
      call heapify_integer(array, perm, n, i)
    end do

    ! Heap sort
    do i = n, 2, -1
      call swap_int(perm(1), perm(i))
      call heapify_integer(array, perm, i - 1, 1)
    end do

  contains

  !> Iterative heapify (non-recursive, pure)
  !| Restore the max-heap property for the subtree rooted at `root`.
  !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
  !| We reorder the permutation vector `perm` (indices) rather than moving
  !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_integer(array, perm, heap_size, root)
      !| Integer input array to sort
      integer(int32), intent(in) :: array(:)
      !| Permutation vector that will be sorted
      integer(int32), intent(inout) :: perm(:)
      !| Size of the heap
      integer(int32), intent(in) :: heap_size
      !| Root index
      integer(int32), intent(in) :: root
      integer(int32) :: current, next_idx, largest_idx

      current = root

      do
        ! Compute left-child index (use largest_idx as left to avoid extra var)
        largest_idx = 2 * current
        ! If there is no left child the subtree is a leaf; nothing to do
        if (largest_idx > heap_size) exit

        ! Compute right-child index (may be out of heap bounds)
        next_idx = largest_idx + 1

        ! Only compare right-child when it actually exists (guarded access)
        if (next_idx <= heap_size) then
          if (array(perm(next_idx)) > array(perm(largest_idx))) then
            ! Right child is larger than left child
            largest_idx = next_idx
          end if
        end if

        ! Compare the selected child with the current node; if the child is
        ! greater, swap permutation indices so the larger value moves up the
        ! heap. We swap entries of `perm` (indices), not the array values.
        if (array(perm(largest_idx)) > array(perm(current))) then
          call swap_int(perm(current), perm(largest_idx))
          current = largest_idx
        else
          exit
        end if
      end do
    end subroutine heapify_integer
  end subroutine heapsort_integer

  !> Internal heapsort implementation for character arrays.
  !| Lexicographic heapsort using string comparison, indirect via `perm`.
  pure subroutine heapsort_character(array, perm)
    !| Character input array to sort
    character(len=*), intent(in)    :: array(:)
    !| Permutation vector that will be sorted
    integer(int32),   intent(inout) :: perm(:)
    integer(int32) :: n, i

    n = size(array)

    ! Build max-heap
    do i = n / 2, 1, -1
      call heapify_character(array, perm, n, i)
    end do

    ! Heap sort
    do i = n, 2, -1
      call swap_int(perm(1), perm(i))
      call heapify_character(array, perm, i - 1, 1)
    end do

  contains

  !> Iterative heapify (non-recursive, pure)
  !| Restore the max-heap property for the subtree rooted at `root`.
  !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
  !| We reorder the permutation vector `perm` (indices) rather than moving
  !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_character(array, perm, heap_size, root)
      !| Character input array to sort
      character(len=*), intent(in) :: array(:)
      !| Permutation vector that will be sorted
      integer(int32),   intent(inout) :: perm(:)
      !| Size of the heap
      integer(int32), intent(in) :: heap_size
      !| Root index
      integer(int32), intent(in) :: root
      integer(int32) :: current, next_idx, largest_idx

      current = root

      do
        ! Compute left-child index (use largest_idx directly as left = 2*current)
        largest_idx = 2 * current
        ! If there is no left child the subtree is a leaf; nothing to do
        if (largest_idx > heap_size) exit

        ! Potential right child index
        next_idx = largest_idx + 1

        ! Only compare right child when it exists to avoid OOB access
        if (next_idx <= heap_size) then
          if (array(perm(next_idx)) > array(perm(largest_idx))) then
            largest_idx = next_idx
          end if
        end if

        ! If the selected child is larger than current, swap permutation
        ! indices so the larger element moves up. We swap entries in `perm`.
        if (array(perm(largest_idx)) > array(perm(current))) then
          call swap_int(perm(current), perm(largest_idx))
          current = largest_idx
        else
          exit
        end if
      end do
    end subroutine heapify_character
  end subroutine heapsort_character

  !> Swap two integer values in-place.
  pure subroutine swap_int(a, b)
    !| First integer to swap
    integer(int32), intent(inout) :: a
    !| Second integer to swap
    integer(int32), intent(inout) :: b
    integer(int32) :: temp
    temp = a; a = b; b = temp
  end subroutine swap_int

  ! Helper: NaN-aware comparisons for real(real64).
  ! We treat NaN as greater than any numeric value so that NaNs end up at the
  ! end of ascending sorts. These helpers define a total-like ordering where
  ! (number) < (NaN) and (NaN) > (number); two NaNs compare equal (neither < nor >).
  !> NaN-aware less-than comparison for real(real64).
!> Finite values are always considered smaller than NaN.
!> Two NaNs are treated as equal (neither less nor greater).
  pure logical function real_less(a, b)
    real(real64), intent(in) :: a
    real(real64), intent(in) :: b

    if (ieee_is_nan(a) .and. ieee_is_nan(b)) then
      real_less = .false.
    else if (ieee_is_nan(a)) then
      real_less = .false.
    else if (ieee_is_nan(b)) then
      real_less = .true.
    else
      real_less = (a < b)
    end if
  end function real_less


  !> NaN-aware greater-than comparison for real(real64).
  !> NaN is treated as greater than any finite number.
  !> Two NaNs are treated as equal.
  pure logical function real_greater(a, b)
    real(real64), intent(in) :: a
    real(real64), intent(in) :: b

    if (ieee_is_nan(a) .and. ieee_is_nan(b)) then
      real_greater = .false.
    else if (ieee_is_nan(a)) then
      real_greater = .true.
    else if (ieee_is_nan(b)) then
      real_greater = .false.
    else
      real_greater = (a > b)
    end if
  end function real_greater

  
  !> Finds the indices of the true values in a logical mask.
  pure subroutine which(mask, n, idx_out, m_max, m_out, ierr)
    !| Logical array of size n.
    logical, intent(in) :: mask(:)
    !| Size of the mask.
    integer(int32), intent(in) :: n
    !| Integer array to store the indices of true values.
    integer(int32), intent(out) :: idx_out(:)
    !| Maximum size of idx_out.
    integer(int32), intent(in) :: m_max
    !| Actual size of idx_out (number of true values found).
    integer(int32), intent(out) :: m_out
    !| Error code: 0=ok, 201=invalid input, 202=empty input
    integer(int32), intent(out) :: ierr

    integer(int32) :: i, count

    ! Initialize error code
    call set_ok(ierr)

    ! Validate inputs
    if (n <= 0 .or. m_max <= 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      m_out = 0
      idx_out = 0
      return
    end if

    if (size(mask) < n .or. size(idx_out) < m_max) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      m_out = 0
      idx_out = 0
      return
    end if

    count = 0
    idx_out = 0  ! Initialize to avoid garbage values
    do i = 1, n
      if (mask(i)) then
        count = count + 1
        if (count <= m_max) then
          idx_out(count) = i
        end if
      end if
    end do
    m_out = count
  end subroutine which

  !> Performs LOESS smoothing on a set of data points.
  !| Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters.
  !| The user must pre-filter data and provide only valid indices in indices_used.
  pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, ierr)
    !| Total number of reference points.
    integer(int32), intent(in) :: n_total
    !| Number of target points to smooth.
    integer(int32), intent(in) :: n_target
    !| Reference x-coordinates.
    real(real64), intent(in) :: x_ref(n_total)
    !| Reference y-coordinates (length n_total).
    real(real64), intent(in) :: y_ref(n_total)
    !| Indices of reference points used for smoothing (only valid indices).
    integer(int32), intent(in) :: indices_used(n_used)
    !| Number of indices actually used for smoothing.
    integer(int32), intent(in) :: n_used
    !| Target x-coordinates to smooth.
    real(real64), intent(in) :: x_query(n_target)
    !| Bandwidth parameter for the kernel.
    real(real64), intent(in) :: kernel_sigma
    !| Cutoff for the kernel.
    real(real64), intent(in) :: kernel_cutoff
    !| Output smoothed values (length n_target).
    real(real64), intent(out) :: y_out(n_target)
    !| Error code: 0=ok, 201=invalid input, 202=empty input
    integer(int32), intent(out) :: ierr

    integer(int32) :: q, i, idx
    real(real64) :: query_x, ref_x, delta, sum_weights, weight
    real(real64) :: min_dist
    integer(int32) :: min_idx
    logical :: exact_match_found, use_kernel

    ! Initialize error code
    call set_ok(ierr)

    ! Input validation
    if (n_total <= 0 .or. n_target <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    if (n_used <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    if (kernel_sigma < 0.0_real64 .or. kernel_cutoff < 0.0_real64) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    ! Validate array sizes
    if (size(x_ref) < n_total .or. size(y_ref) < n_total .or. &
        size(indices_used) < n_used .or. size(x_query) < n_target .or. &
        size(y_out) < n_target) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    ! Validate indices are within bounds
    do i = 1, n_used
      if (indices_used(i) < 1 .or. indices_used(i) > n_total) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        y_out = 0.0_real64
        return
      end if
    end do

    ! Check if we should use kernel smoothing
    use_kernel = (kernel_sigma > 0.0_real64)

    do q = 1, n_target
      query_x = x_query(q)
      sum_weights = 0.0_real64
      y_out(q) = 0.0_real64
      min_dist = huge(1.0_real64)
      min_idx = indices_used(1)
      exact_match_found = .false.

      ! Process all reference points
      do i = 1, n_used
        idx = indices_used(i)
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        
        ! Check for exact match
        if (delta == 0.0_real64) then
          y_out(q) = y_ref(idx)
          exact_match_found = .true.
          exit
        end if
        
        ! Track closest point for potential fallback
        if (delta < min_dist) then
          min_dist = delta
          min_idx = idx
        end if
        
        ! Apply kernel smoothing if enabled and within cutoff
        if (use_kernel .and. delta <= kernel_cutoff * kernel_sigma) then
          weight = exp(-(delta / kernel_sigma)**2)
          sum_weights = sum_weights + weight
          y_out(q) = y_out(q) + weight * y_ref(idx)
        end if
      end do

      ! Finalize result if no exact match was found
      if (.not. exact_match_found) then
        if (sum_weights > 0.0_real64) then
          ! We have weighted average from kernel smoothing
          y_out(q) = y_out(q) / sum_weights
        else
          ! Fallback: use nearest neighbor
          y_out(q) = y_ref(min_idx)
        end if
      end if
    end do
  end subroutine loess_smooth_2d

  !> Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
  !| Returns the sorted unique values and their cumulative frequencies in [0,1].
  !| Assumes perm is already sorted by values[perm]. Caller controls sorting algorithm.
  !| The number of unique values can be determined by finding the last non-zero cdf_value.
  subroutine compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
    !| Array of observed data values (e.g., contributions or spikes).
    real(real64), intent(in) :: values(n_values)
    !| Number of values in the input array.
    integer(int32), intent(in) :: n_values
    !| Pre-sorted permutation indices (must be sorted by values[perm]).
    integer(int32), intent(in) :: perm(n_values)
    !| Sorted unique data values.
    real(real64), intent(out) :: unique_values(n_values)
    !| Corresponding cumulative frequencies between 0 and 1.
    real(real64), intent(out) :: cdf_values(n_values)
    !| Number of unique values found (actual size of output arrays)
    integer(int32), intent(out) :: n_unique
    !| Error code: 0=ok, 201=invalid input, 202=empty input
    integer(int32), intent(out) :: ierr

    integer(int32) :: i
    real(real64) :: current_val, cumulative_count

    ! Initialize error code and outputs
    call set_ok(ierr)
    unique_values = 0.0_real64
    cdf_values = 0.0_real64

    ! Input validation
    if (n_values <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if

    if (size(values) < n_values .or. size(perm) < n_values) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! Identify unique values and compute cumulative frequencies
    n_unique = 0
    cumulative_count = 0.0_real64
    current_val = -huge(1.0_real64)  ! Initialize to lowest possible value

    do i = 1, n_values
      ! Check if this is a new unique value (exact comparison)
      if (values(perm(i)) /= current_val) then
        ! New unique value found
        current_val = values(perm(i))
        n_unique = n_unique + 1
        
        if (n_unique > size(unique_values)) then
          ! Output array is too small
          call set_err_once(ierr, ERR_INVALID_INPUT)
          return
        end if
        
        unique_values(n_unique) = current_val
      end if
      
      ! Update cumulative count and set CDF value
      cumulative_count = cumulative_count + 1.0_real64
      cdf_values(n_unique) = cumulative_count / real(n_values, real64)
    end do
  end subroutine compute_edf

  !> Helper routine that sorts and calls compute_edf.
  !| Allocates workspace internally and performs sorting before computing EDF.
  !| Use this for convenience; use compute_edf directly for custom sorting.
  subroutine compute_edf_alloc(values, n_values, unique_values, cdf_values, n_unique, ierr)
    !| Array of observed data values (e.g., contributions or spikes).
    real(real64), intent(in) :: values(n_values)
    !| Number of values in the input array.
    integer(int32), intent(in) :: n_values
    !| Sorted unique data values.
    real(real64), intent(out) :: unique_values(n_values)
    !| Corresponding cumulative frequencies between 0 and 1.
    real(real64), intent(out) :: cdf_values(n_values)
    !| Number of unique values found (actual size of output arrays)
    integer(int32), intent(out) :: n_unique
    !| Error code: 0=ok, 201=invalid input, 202=empty input
    integer(int32), intent(out) :: ierr

    ! Local workspace arrays with explicit size
    integer(int32) :: perm(n_values)
    integer(int32) :: stack_left(n_values)
    integer(int32) :: stack_right(n_values)
    integer(int32) :: i

    ! Initialize permutation vector to [1, 2, 3, ..., n_values]
    do i = 1, n_values
      perm(i) = i
    end do

    ! Sort values using the permutation vector
    call sort_array(values, perm, stack_left, stack_right)

    ! Compute EDF with sorted permutation
    call compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
  end subroutine compute_edf_alloc

  !> Calculate the percentile of an array given a sorted permutation.
  !! Uses linear interpolation between adjacent values.
  pure subroutine calc_percentile(array, permutation, percentile, value, ierr)
    real(real64), intent(in) :: array(:)
    !! input array
    integer(int32), intent(in) :: permutation(:)
    !! permutation vector representing sorted order
    real(real64), intent(in) :: percentile
    !! desired percentile (0-100)
    real(real64), intent(out) :: value
    !! output percentile value
    integer(int32), intent(out) :: ierr
    !! Error code
    
    integer(int32) :: n, lower_index
    real(real64) :: index, fraction, lower_value, upper_value
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Input validation
    n = size(array)
    if (n == 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      value = 0.0_real64
      return
    end if
    
    if (size(permutation) /= n) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      value = 0.0_real64
      return
    end if
    
    if (percentile < 0.0_real64 .or. percentile > 100.0_real64) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      value = 0.0_real64
      return
    end if
    
    ! Handle single element case
    if (n == 1) then
      value = array(1)
      return
    end if
    
    ! Calculate the fractional index using linear interpolation method
    ! This follows the method used in numpy.percentile with interpolation='linear'
    index = (percentile / 100.0_real64) * real(n - 1, real64) + 1.0_real64
    
    lower_index = floor(index)
    fraction = index - real(lower_index, real64)
    
    ! Handle edge cases for indices
    if (lower_index < 1) then
      value = array(permutation(1))  ! Smallest value in sorted order
    else if (lower_index >= n) then
      value = array(permutation(n))  ! Largest value in sorted order
    else
      ! Linear interpolation between adjacent values using permuted indices
      lower_value = array(permutation(lower_index))
      upper_value = array(permutation(lower_index + 1))
      value = lower_value + fraction * (upper_value - lower_value)
    end if
  end subroutine calc_percentile

  !> Calculate the percentile of an array, allocating necessary arrays when no sorting permutation is given
  !! @note This subroutine uses quicksort internally which may cause a spike in memory usage for large arrays.
  subroutine calc_percentile_alloc(array, percentile, value, ierr)
    use tox_errors, only: ERR_EMPTY_INPUT, ERR_ALLOC_FAIL, set_ok, set_err, is_err
    real(real64), intent(in) :: array(:)
    !! Input array
    real(real64), intent(in) :: percentile
    !! Desired percentile (0-100)
    real(real64), intent(out) :: value
    !! Output percentile value
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: n, i
    integer(int32), allocatable :: perm(:)
    integer(int32), allocatable :: stack_left(:), stack_right(:)

    n = size(array)
    ! Initialize error
    call set_ok(ierr)

    if (n == 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      value = 0.0_real64
      return
    end if

    ! Allocate permutation and stacks
    allocate(perm(n), stack_left(n), stack_right(n), stat=ierr)
    if (is_err(ierr)) then
      call set_err(ierr, ERR_ALLOC_FAIL)
      value = 0.0_real64
      return
    end if

    ! Initialize permutation to identity
    perm = [(i, i=1,n)]

    ! Sort the array indirectly
    call sort_real(array, perm, stack_left, stack_right)

    ! Calculate percentile using sorted permutation
    call calc_percentile(array, perm, percentile, value, ierr)
  end subroutine calc_percentile_alloc

  pure function get_median(dist, perm, k_real) result(res)
    real(real64), intent(in) :: dist(:)
    integer(int32), intent(in) :: perm(:), k_real
    real(real64) :: res
    ! El índice 1 es self, la mediana de los vecinos empieza después
    res = dist(perm(1 + (k_real+1)/2))
  end function

  pure function count_valid(dist, perm, k_total) result(res)
    real(real64), intent(in) :: dist(:)
    integer(int32), intent(in) :: perm(:), k_total
    integer(int32) :: res, j
    res = 0
    do j = 2, k_total
      if (dist(perm(j)) < 1.0e20_real64) res = res + 1
    end do
  end function

  pure function calculate_mean(arr, n) result(res)
    real(real64), intent(in) :: arr(:)
    integer(int32), intent(in) :: n
    real(real64) :: res
    res = sum(arr) / real(n, real64)
  end function

  !> This subroutine calculates the roughness of a smoothed field based on the
  !! geometric neighborhood structure. Roughness quantifies the average weighted
  !! disagreement between neighboring points, providing a measure of local oscillations
  !! or irregularities in the field.
  pure subroutine compute_roughness(y, sigma, neighbors, distances, &
                             n_points, n_dims, k_neighbors, roughness)
      integer(int32), intent(in) :: n_points       
      !| Number of points in the dataset
      integer(int32), intent(in) :: n_dims         
      !| Number of dimensions for each point
      integer(int32), intent(in) :: k_neighbors    
      !| Number of nearest neighbors per point
      real(real64), intent(in) :: y(n_dims, n_points)       
      !| Smoothed field values for each point
      real(real64), intent(in) :: sigma(n_points)           
      !| Adaptive bandwidth for each point
      integer(int32), intent(in) :: neighbors(k_neighbors, n_points)  
      !| Indices of k-nearest neighbors
      real(real64), intent(in) :: distances(k_neighbors, n_points)  
      !| Distances to k-nearest neighbors
      real(real64), intent(out) :: roughness     
      !| Computed roughness value

      ! --- Local variables ---
      integer(int32) :: i, j_idx, neighbor_id, m   
      !| Loop indices and neighbor identifiers
      real(real64) :: d2, w_ij, w_ji, W_sym      
      !| Distance squared and weights
      real(real64) :: sum_diff2           
      !| Squared difference of smoothed field
      real(real64) :: R_acc, Z_acc               
      !| Accumulators for roughness and normalization
      real(real64) :: inv_2s2(n_points)
      !| Inverse of 2 * sigma^2 for points i and j

      ! Precompute the inverse factor for sigma of point i
      do i = 1_int32, n_points
          if (sigma(i) > 1.0e-12_real64) then
              inv_2s2(i) = 1.0_real64 / (2.0_real64 * sigma(i)**2)
          else
              inv_2s2(i) = 0.0_real64
          end if
      end do

      ! 1. Initialize accumulators
      R_acc = 0.0_real64  ! Total weighted local variation
      Z_acc = 0.0_real64  ! Total weight mass for normalization

      ! 2. Iterate over each point i
      do i = 1_int32, n_points
          
          if (inv_2s2(i) <= 0.0_real64) cycle

          ! Iterate over each neighbor j of point i
          do j_idx = 1_int32, k_neighbors
              neighbor_id = neighbors(j_idx, i)  ! Neighbor index

              if (neighbor_id < 1_int32 .or. neighbor_id > n_points) cycle

              ! Skip self-loops if the kNN includes the point itself
              if (neighbor_id == i) cycle

              if (inv_2s2(neighbor_id) <= 0.0_real64) cycle

              ! (a) Compute directed weights
              ! Use the precomputed distance d_ij
              d2 = distances(j_idx, i)**2

              w_ij = exp(-d2 * inv_2s2(i))
              w_ji = exp(-d2 * inv_2s2(neighbor_id))

              ! (b) Compute symmetric weight (geometric mean)
              W_sym = sqrt(w_ij * w_ji)

              ! (c) Compute squared difference of the smoothed field (diff2)
              sum_diff2 = 0.0_real64
              do m = 1_int32, n_dims
                  sum_diff2 = sum_diff2 + (y(m, i) - y(m, neighbor_id))**2
              end do

              ! (d) Accumulate contributions
              R_acc = R_acc + (W_sym * sum_diff2)  ! Weighted local variation
              Z_acc = Z_acc + W_sym                ! Total weight mass
          end do
      end do

      ! 3. Final normalization
      if (Z_acc > 1.0e-15) then
          roughness = R_acc / Z_acc  ! Normalize total variation by total weight
      else
          roughness = 0.0_real64  ! Default to zero if no weights are present
      end if

  end subroutine compute_roughness

  !> This subroutine calculates the Root Mean Square Error (RMSE) between
  !! the current smoothed field (y) and the original input field (y0).
  !! It serves as a measure of global distortion or fidelity loss.
  pure subroutine compute_rmse(y0, y, n_points, n_dims, rmse)      
      integer(int32), intent(in) :: n_points
      !| Number of points 
      integer(int32), intent(in) :: n_dims
      !| Number of dimensions
      real(real64), intent(in) :: y0(n_dims, n_points) 
      !| Original data
      real(real64), intent(in) :: y(n_dims, n_points)  
      !| Current smoothed field
      real(real64), intent(out) :: rmse                
      !| Output scalar

      ! --- Local variables ---
      integer(int32) :: i, m
      real(real64) :: s_acc, diff_val

      s_acc = 0.0_real64

      ! i as outer loop + m as inner loop = Stride-1 access (Fastest in Fortran)
      do i = 1_int32, n_points
          do m = 1_int32, n_dims
              diff_val = y(m, i) - y0(m, i)
              s_acc = s_acc + (diff_val * diff_val) ! Faster than **2
          end do
      end do

      if (n_points > 0_int32) then
          ! Normalizing by n_points as per user specification (Euclidean drift)
          rmse = sqrt(s_acc / real(n_points, real64))
      else
          rmse = 0.0_real64
      end if

  end subroutine compute_rmse
  
  !> Calculates the Coverage metric and its penalty using a precomputed reference diameter
  !> and frozen anchor indices (i1,i2) computed at t0.
  pure subroutine compute_coverage(d0, i1, i2, y, n_points, n_dims, coverage, coverage_penalty)
      real(real64), intent(in) :: d0
      !| Precomputed original diameter (from t0)
      integer(int32), intent(in) :: i1, i2
      !| Frozen anchor indices (from t0)
      integer(int32), intent(in) :: n_points
      !| Number of points
      integer(int32), intent(in) :: n_dims
      !| Number of dimensions
      real(real64), intent(in) :: y(n_dims, n_points)
      !| Current field
      real(real64), intent(out) :: coverage
      !| Proportion of points covered by the algorithm. (higher is better)
      real(real64), intent(out) :: coverage_penalty
      !| Associated with the coverage, used to adjust or regularize the result (lower is better)

      ! --- Local variables ---
      real(real64) :: dt_sq, dt, diff
      integer(int32) :: m
      real(real64), parameter :: eps = 1.0e-12_real64

      ! 1. Compute current anchor distance (frozen anchors)
      dt_sq = 0.0_real64
      do m = 1, n_dims
          diff = y(m, i1) - y(m, i2)
          dt_sq = dt_sq + diff * diff
      end do
      dt = sqrt(max(0.0_real64, dt_sq))

      ! 2. Compute Coverage ratio using the passed d0
      if (d0 > eps) then
          coverage = dt / d0
      else
          ! Degenerate cases: if original was a point
          if (dt < eps) then
              coverage = 1.0_real64
          else
              coverage = 0.0_real64
          end if
      end if

      ! 3. Compute Penalty (lower is better)
      coverage_penalty = 1.0_real64 / (coverage + eps)

  end subroutine compute_coverage


  !> 2-sweep farthest-point heuristic for approximate diameter.
  !> Also returns the anchor indices (i1,i2) defining the approximate diameter.
  pure subroutine get_approx_diameter(field, n, dim, diameter, i1_out, i2_out)
      integer(int32), intent(in) :: n
      !| Number of points in the field.
      integer(int32), intent(in) :: dim
      !| Dimensionality of the field
      real(real64), intent(in) :: field(dim, n)
      !| Input array representing the field
      real(real64), intent(out) :: diameter
      !| diameter of the field
      integer(int32), intent(out) :: i1_out, i2_out
      !| Anchor indices defining the approximate diameter

      integer(int32) :: j, m
      real(real64) :: max_dist_sq, current_dist_sq, diff

      if (n < 2_int32) then
          diameter = 0.0_real64
          i1_out = 1
          i2_out = 1
          return
      end if

      ! Step 1: Find i1 farthest from index 1 (O(n))
      i1_out = 1
      max_dist_sq = -1.0_real64
      do j = 1, n
          current_dist_sq = 0.0_real64
          do m = 1, dim
              diff = field(m, j) - field(m, 1)
              current_dist_sq = current_dist_sq + diff * diff
          end do
          if (current_dist_sq > max_dist_sq) then
              max_dist_sq = current_dist_sq
              i1_out = j
          end if
      end do

      ! Step 2: Find i2 farthest from i1 (O(n))
      i2_out = i1_out
      max_dist_sq = -1.0_real64
      do j = 1, n
          current_dist_sq = 0.0_real64
          do m = 1, dim
              diff = field(m, j) - field(m, i1_out)
              current_dist_sq = current_dist_sq + diff * diff
          end do
          if (current_dist_sq > max_dist_sq) then
              max_dist_sq = current_dist_sq
              i2_out = j
          end if
      end do

      if (max_dist_sq > 0.0_real64) then
          diameter = sqrt(max_dist_sq)
      else
          diameter = 0.0_real64
      end if
  end subroutine get_approx_diameter

  !> Computes the composite smoothing score (St) using a geometric mean.
  !! Balances roughness, fidelity (RMSE), and coverage.
  !! Lower is better.
  pure subroutine compute_smoothing_score(roughness_t, roughness_0, &
                                          rmse_t, diameter_0, &
                                          coverage_penalty_t, epsilon, &
                                          score_t)
    
      real(real64), intent(in) :: roughness_t
      !| Current roughness at iteration t.
      real(real64), intent(in) :: roughness_0
      !| Reference roughness of the original noisy data (t=0). Used for normalization.
      real(real64), intent(in) :: rmse_t
      !| Root Mean Square Error between current field and original input.
      real(real64), intent(in) :: diameter_0
      !| Approximate geometric diameter of original data. Acts as the global length scale.
      real(real64), intent(in) :: coverage_penalty_t
      !| Dimensionless penalty (1/coverage). Tracks global contraction or collapse.
      real(real64), intent(in) :: epsilon
      !| Safety constant to prevent division by zero in normalizations.
      real(real64), intent(out) :: score_t
      !| Composite dimensionless score (geometric mean). Lower values indicate better trade-offs.

      ! --- Local variables ---
      real(real64) :: r_norm, e_norm, combined_product
      real(real64) :: Rt, R0, Et, D0, pt, eps
      real(real64), parameter :: one_third = 0.3333333333333333_real64

      ! Defensive clamps (these quantities should be non-negative)
      Rt  = max(0.0_real64, roughness_t)
      R0  = max(0.0_real64, roughness_0)
      Et  = max(0.0_real64, rmse_t)
      D0  = max(0.0_real64, diameter_0)
      pt  = max(0.0_real64, coverage_penalty_t)
      eps = max(0.0_real64, epsilon)

      ! 1. Normalize Roughness (starts at ~1.0)
      r_norm = Rt / (R0 + eps)

      ! 2. Normalize RMSE (starts at 0.0)
      ! Normalized by diameter to make it dimensionless and scale-independent
      e_norm = Et / (D0 + eps)

      ! 3. Composite product
      ! We use (1 + e_norm) to prevent the product from being zero at t=0
      combined_product = r_norm * (1.0_real64 + e_norm) * pt

      ! 4. Geometric Mean (Cube root)
      ! Lowering any factor lowers the score.
      if (combined_product > 0.0_real64) then
          score_t = combined_product**one_third
      else
          ! If something degenerates (e.g., pt=0), return a large score rather than "perfect"
          score_t = huge(1.0_real64)
      end if

  end subroutine compute_smoothing_score



end module f42_utils

! === R WRAPPERS ===

!> R wrapper for loess_smooth_2d.
!| Direct wrapper - user must pre-filter indices in R before calling.
subroutine loess_smooth_2d_r(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr)
  use f42_utils, only: loess_smooth_2d
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  !| Total number of reference points.
  integer(int32), intent(in) :: n_total
  !| Number of target points to smooth.
  integer(int32), intent(in) :: n_target
  !| Reference x-coordinates.
  real(real64), intent(in) :: x_ref(n_total)
  !| Reference y-coordinates (length n_total).
  real(real64), intent(in) :: y_ref(n_total)
  !| Indices of reference points used for smoothing (pre-filtered).
  integer(int32), intent(in) :: indices_used(n_used)
  !| Number of indices actually used for smoothing.
  integer(int32), intent(in) :: n_used
  !| Target x-coordinates to smooth.
  real(real64), intent(in) :: x_query(n_target)
  !| Bandwidth parameter for the kernel.
  real(real64), intent(in) :: kernel_sigma
  !| Cutoff for the kernel.
  real(real64), intent(in) :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(real64), intent(out) :: y_out(n_target)
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(int32), intent(out) :: ierr
  
  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr)
end subroutine loess_smooth_2d_r

! === C WRAPPERS ===

!> C wrapper for which.
!| Converts integer mask to logical and calls which.
subroutine which_c(mask, n, idx_out, m_max, m_out, ierr) bind(C, name="which_c")
  use, intrinsic :: iso_c_binding, only: c_int
  use, intrinsic :: iso_fortran_env, only: int32
  use f42_utils, only: which
  implicit none
  !| Size of the mask.
  integer(c_int), intent(in), value :: n
  !| Maximum size of idx_out.
  integer(c_int), intent(in), value :: m_max
  !| Integer mask array (0/1 values).
  integer(c_int), intent(in) :: mask(n)
  !| Output array for indices of true values.
  integer(c_int), intent(out) :: idx_out(m_max)
  !| Actual size of idx_out (number of true values found).
  integer(c_int), intent(out) :: m_out
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr
  logical :: mask_f(n)
  integer(int32) :: i, ierr_f
  do i = 1, n
    mask_f(i) = (mask(i) /= 0)
  end do
  call which(mask_f, n, idx_out, m_max, m_out, ierr_f)
  ierr = ierr_f
end subroutine which_c

!> C wrapper for loess_smooth_2d.
!| Direct wrapper - user must pre-filter indices in C before calling.
subroutine loess_smooth_2d_c(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr) bind(C, name="loess_smooth_2d_c")
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use, intrinsic :: iso_fortran_env, only: int32
  use f42_utils, only: loess_smooth_2d
  implicit none
  !| Total number of reference points.
  integer(c_int), intent(in), value :: n_total
  !| Number of target points to smooth.
  integer(c_int), intent(in), value :: n_target
  !| Reference x-coordinates.
  real(c_double), intent(in) :: x_ref(n_total)
  !| Reference y-coordinates (length n_total).
  real(c_double), intent(in) :: y_ref(n_total)
  !| Indices of reference points used for smoothing (pre-filtered).
  integer(c_int), intent(in) :: indices_used(n_used)
  !| Number of indices actually used for smoothing.
  integer(c_int), intent(in), value :: n_used
  !| Target x-coordinates to smooth.
  real(c_double), intent(in) :: x_query(n_target)
  !| Bandwidth parameter for the kernel.
  real(c_double), intent(in), value :: kernel_sigma
  !| Cutoff for the kernel.
  real(c_double), intent(in), value :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(c_double), intent(out) :: y_out(n_target)
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr

  integer(int32) :: ierr_f
  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr_f)
  ierr = ierr_f

end subroutine loess_smooth_2d_c

!> IMPORTANT - When using these C wrapper functions, no copies of the arrays will be created.
!| The Fortran routine will operate directly on the memory provided by the caller.

!> C wrapper for compute_edf.
!| Allocates workspace internally and exposes a simple interface with C types.
!| The number of unique values is returned via n_unique output parameter.
subroutine compute_edf_c(values, n_values, unique_values, cdf_values, n_unique, ierr) &
  bind(C, name="compute_edf_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_utils, only: compute_edf_alloc
  implicit none
  !| Number of values in the input array.
  integer(c_int), intent(in), value :: n_values
  !| Array of observed data values (e.g., contributions or spikes).
  real(c_double), intent(in) :: values(n_values)
  !| Sorted unique data values (sized to n_values).
  real(c_double), intent(out) :: unique_values(n_values)
  !| Corresponding cumulative frequencies between 0 and 1 (sized to n_values).
  real(c_double), intent(out) :: cdf_values(n_values)
  !| Number of unique values found.
  integer(c_int), intent(out) :: n_unique
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr

  call compute_edf_alloc(values, n_values, unique_values, cdf_values, n_unique, ierr)
end subroutine compute_edf_c

!> Expert C wrapper for compute_edf.
!| Direct interface to compute_edf for users who have already sorted their data
!| or have a custom permutation vector. This skips the internal sorting step
!| for better performance when the caller has full control.
subroutine compute_edf_expert_c(values, n_values, perm, unique_values, cdf_values, n_unique, ierr) &
    bind(C, name="compute_edf_expert_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_utils, only: compute_edf
  implicit none
  !| Number of values in the input array.
  integer(c_int), intent(in), value :: n_values
  !| Array of observed data values (e.g., contributions or spikes).
  real(c_double), intent(in) :: values(n_values)
  !| Pre-sorted permutation indices (must be sorted by values[perm]).
  !| Caller is responsible for sorting this array before calling.
  integer(c_int), intent(in) :: perm(n_values)
  !| Sorted unique data values (sized to n_values).
  real(c_double), intent(out) :: unique_values(n_values)
  !| Corresponding cumulative frequencies between 0 and 1 (sized to n_values).
  real(c_double), intent(out) :: cdf_values(n_values)
  !| Number of unique values found.
  integer(c_int), intent(out) :: n_unique
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr

  call compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
end subroutine compute_edf_expert_c
