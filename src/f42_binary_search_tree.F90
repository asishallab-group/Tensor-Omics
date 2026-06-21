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
  
  !> category: C-interface
  !| Build the BST index by sorting indices using values in x.
  pure subroutine build_bst_index(values, num_values, sorted_indices, tmp_left_stack, tmp_right_stack, ierr)
    integer(int32), intent(in) :: num_values           
    !! Number of elements in values array
    real(real64), intent(in) :: values(num_values)      
    !! Input real array to be indexed
    integer(int32), intent(out) :: sorted_indices(num_values)  
    !! Output permutation index
    integer(int32), intent(out) :: tmp_left_stack(num_values)    
    !! Manual stack for left indices
    integer(int32), intent(out) :: tmp_right_stack(num_values)   
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
    call sort_array(values, sorted_indices, tmp_left_stack, tmp_right_stack)
  end subroutine build_bst_index

  !> category: C-interface
  !| Get the value at the sorted position.
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

  !> category: C-interface
  !| Perform a 1D range query over the sorted index.
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
    !! DM_RESULT_SIZE_IS(num_matches)
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
