!> Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.

module f42_utils
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once
  implicit none

  public :: sort_real, sort_integer, sort_character
  public :: sort_array

  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array
contains

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
    call quicksort_real(array, perm, size(array), stack_left, stack_right)
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
    call quicksort_int(array, perm, size(array), stack_left, stack_right)
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
    call quicksort_char(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_character

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

  !> Swap two integer values in-place.
  pure subroutine swap_int(a, b)
    !| First integer to swap
    integer(int32), intent(inout) :: a
    !| Second integer to swap
    integer(int32), intent(inout) :: b
    integer(int32) :: temp
    temp = a; a = b; b = temp
  end subroutine swap_int

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
  use iso_c_binding
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
  use iso_c_binding, only : c_int, c_double
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
