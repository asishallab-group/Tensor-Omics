#include "macros.h"

!! Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!!
!! This module provides routines to build a BST index (via sorting), access sorted values,
!! and perform efficient range queries over a real-valued array.
module f42_binary_search_tree
  use safeguard
  use f42_utils, only: sort_array
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_DIM_MISMATCH, is_ok, set_err_once, set_ok, validate_dimension_size
  implicit none
  public :: build_bst_index, get_sorted_value, bst_range_query
contains
  
  !> Build the BST index by sorting indices using values in x.
  pure subroutine build_bst_index(values, num_values, sorted_indices, left_stack, right_stack, ierr)
    integer(int32), intent(in) :: num_values           
    !! Number of elements in values array
    real(real64), intent(in) :: values(num_values)      
    !! Input real array to be indexed
    integer(int32), intent(out) :: sorted_indices(num_values)  
    !! Output permutation index
    integer(int32), intent(out) :: left_stack(num_values)    
    !! Manual stack for left indices
    integer(int32), intent(out) :: right_stack(num_values)   
    !! Manual stack for right indices
    integer(int32), intent(out) :: ierr                   
    !! Error code
    integer(int32) :: idx

    call set_ok(ierr)
    
    call validate_dimension_size(num_values, ierr)
    if(.not. is_ok(ierr)) return

    do idx = 1, num_values
      sorted_indices(idx) = idx
    end do
    call sort_array(values, sorted_indices, left_stack, right_stack)
  end subroutine build_bst_index

  !> Get the value at the sorted position.
  function get_sorted_value(values, sorted_indices, position, ierr) result(sorted_value)
    real(real64), intent(in) :: values(:)              
    !! Input real array
    integer(int32), intent(in) :: sorted_indices(:)    
    !! Permutation index array
    integer(int32), intent(in) :: position             
    !! Sorted position (1-based)
    integer(int32), intent(out) :: ierr                
    !! Error code
    real(real64) :: sorted_value

    call set_ok(ierr)
    
    ! Input validation
    if (position < 1 .or. position > size(sorted_indices)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    if (sorted_indices(position) < 1 .or. sorted_indices(position) > size(values)) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    sorted_value = values(sorted_indices(position))
  end function get_sorted_value

  !> Perform a 1D range query over the sorted index.
  pure subroutine bst_range_query(values, sorted_indices, num_values, lower_bound, upper_bound, &
                            output_indices, num_matches, ierr)

    integer(int32), intent(in) :: num_values           
    !! Number of elements                        
    real(real64), intent(in) :: values(num_values)      
    !! Input real array
    integer(int32), intent(in) :: sorted_indices(num_values)  
    !! Permutation index array (sorted) 
    real(real64), intent(in) :: lower_bound            
    !! Lower bound of range (inclusive)
    real(real64), intent(in) :: upper_bound            
    !! Upper bound of range (inclusive)
    integer(int32), intent(out) :: output_indices(num_values)  
    !! Output array of matching indices
    integer(int32), intent(out) :: num_matches         
    !! Number of matches found
    integer(int32), intent(out) :: ierr                
    !! Error code
    integer(int32) :: idx

    call set_ok(ierr)
    
    call validate_dimension_size(num_values, ierr)
    if(.not. is_ok(ierr)) return
    
    if (lower_bound > upper_bound) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    num_matches = 0
    do idx = 1, num_values
      if (values(sorted_indices(idx)) >= lower_bound .and. &
          values(sorted_indices(idx)) <= upper_bound) then
        num_matches = num_matches + 1
        output_indices(num_matches) = sorted_indices(idx)
      else if (values(sorted_indices(idx)) > upper_bound) then
        exit
      end if
    end do
  end subroutine bst_range_query

end module f42_binary_search_tree



!> Wrapper using C for getting range query usable by python
pure subroutine bst_range_query_C(values, sorted_indices, num_values, lower_bound, upper_bound, &
                            output_indices, num_matches, ierr) bind(C, name='bst_range_query_C')
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use f42_binary_search_tree, only: bst_range_query
  M_USE_NULL_VALIDATION
  implicit none
  real(c_double), intent(in), target :: values(num_values)    
  !! Input real array (C-style)
  integer(c_int), intent(in), target :: sorted_indices(num_values)  
  !! Permutation index array (C-style)
  integer(c_int), intent(in), target :: num_values               
  !! Number of elements
  real(c_double), intent(in), target :: lower_bound              
  !! Lower bound of range
  real(c_double), intent(in), target :: upper_bound              
  !! Upper bound of range
  integer(c_int), intent(out), target :: output_indices(num_values)  
  !! Output array (C-style)
  integer(c_int), intent(out), target :: num_matches        
  !! Number of matches found
  integer(c_int), intent(out), target :: ierr                
  !! Error code

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(num_values)
  M_CHECK_NON_NULL(lower_bound)
  M_CHECK_NON_NULL(upper_bound)
  M_CHECK_NON_NULL(values)
  M_CHECK_NON_NULL(sorted_indices)
  M_CHECK_NON_NULL(output_indices)
  M_CHECK_NON_NULL(num_matches)

  call bst_range_query(values, sorted_indices, num_values, lower_bound, upper_bound, &
                      output_indices, num_matches, ierr)
end subroutine bst_range_query_C

!> Wrapper using C for building BST index usable by python
pure subroutine build_bst_index_C(values, num_values, sorted_indices, left_stack, right_stack, ierr) &
                            bind(C, name='build_bst_index_C')
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use f42_binary_search_tree
  M_USE_NULL_VALIDATION
  implicit none
  integer(c_int), intent(in), target :: num_values               
  !! Number of elements
  real(c_double), intent(in), target :: values(num_values)    
  !! Input real array (C-style)
  integer(c_int), intent(out), target :: sorted_indices(num_values)  
  !! Output permutation index (C-style)
  integer(c_int), intent(out), target :: left_stack(num_values)    
  !! Manual stack for left indices
  integer(c_int), intent(out), target :: right_stack(num_values)   
  !! Manual stack for right indices
  integer(c_int), intent(out), target :: ierr                
  !! Error code

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(num_values)
  M_CHECK_NON_NULL(values)
  M_CHECK_NON_NULL(sorted_indices)
  M_CHECK_NON_NULL(left_stack)
  M_CHECK_NON_NULL(right_stack)

  call build_bst_index(values, num_values, sorted_indices, left_stack, right_stack, ierr)
end subroutine build_bst_index_C