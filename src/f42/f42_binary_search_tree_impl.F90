#include <src/macros.h>

!> Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!| This module provides routines to build a BST index (via sorting), access sorted values,
!| and perform range queries over a real-valued array using the sorted index.
module f42_binary_search_tree_impl
    use f42_sort_impl, only: sort_array_heapsort, init_perm, binary_search_insertion
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err, set_ok, validate_in_range_int
    M_IMPLICIT_NONE
    public :: build_bst_index_impl, get_sorted_value, bst_range_query_impl
contains

    !> summary: Build the BST index by sorting indices using values in x
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine build_bst_index_impl(values, n_values, sorted_indices)
        integer(int32), intent(in) :: n_values
            !! Number of elements in values array
        real(real64), intent(in) :: values(n_values)
            !! Input real array to be indexed
        integer(int32), intent(out) :: sorted_indices(n_values)
            !! Output permutation index

        call init_perm(sorted_indices)
        call sort_array_heapsort(values, sorted_indices)
    end subroutine build_bst_index_impl

    !> AUTHOR_AARON_SCHROEDER
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
        call validate_in_range_int(position, ierr, min=1_int32, max=size(sorted_indices, kind=int32), arg_pos=3_int32)
        if (is_err(ierr)) return
        call validate_in_range_int(sorted_indices(position), ierr, min=1_int32, max=size(values, kind=int32), arg_pos=2_int32)
        if (is_err(ierr)) return

        sorted_value = values(sorted_indices(position))
    end function get_sorted_value

    !> summary: Perform a 1D range query over the sorted index
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine bst_range_query_impl(values, values_perm, n_values, lower_bound, upper_bound, &
                                         output_indices, n_matches)

        integer(int32), intent(in) :: n_values
            !! Number of elements
        real(real64), intent(in) :: values(n_values)
            !! Input real array
        integer(int32), intent(in) :: values_perm(n_values)
            !! Permutation of `values` in ascending order -- the BST index. The allocating entry
            !! point builds and heapsorts it for you; the expert one takes whatever order you
            !! supply, so a caller that already holds one from
            !! [[f42_binary_search_tree(module):build_bst_index]] can reuse it across queries.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_values)
        real(real64), intent(in) :: lower_bound
            !! Lower bound of range (inclusive)
            !! DM_MAX(upper_bound)
        real(real64), intent(in) :: upper_bound
            !! Upper bound of range (inclusive)
        integer(int32), intent(out) :: output_indices(n_values)
            !! Output array of matching indices.
            !! DM_RESULT_SIZE_IS(n_matches)
        integer(int32), intent(out) :: n_matches
            !! Number of matches found
        integer(int32) :: i_value

        n_matches = 0

        ! find index of first value in range
        i_value = binary_search_insertion(values, values_perm, lower_bound)
        do while (i_value <= n_values)
            ! check if value is in range and add to output, else exit loop
            if (values(values_perm(i_value)) <= upper_bound) then
                n_matches = n_matches + 1
                output_indices(n_matches) = values_perm(i_value)
                i_value = i_value + 1
            else
                exit
            end if
        end do
    end subroutine bst_range_query_impl

end module f42_binary_search_tree_impl
