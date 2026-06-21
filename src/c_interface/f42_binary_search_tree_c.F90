#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[f42_binary_search_tree(module)]]
module f42_binary_search_tree_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[f42_binary_search_tree(module):get_sorted_value(function)]]
    !| Get the value at the sorted position.
    subroutine get_sorted_value_c(&
            values,&
            n_values_elements,&
            sorted_indices,&
            n_sorted_indices_elements,&
            position,&
            ierr,&
            sorted_value&
            ) bind(C, name="get_sorted_value_c")
        use f42_binary_search_tree, only: get_sorted_value
        use f42_binary_search_tree
        integer(c_int), intent(in), target :: n_values_elements
            !!  Size of the 1. dimension/extent of `values`
        integer(c_int), intent(in), target :: n_sorted_indices_elements
            !!  Size of the 1. dimension/extent of `sorted_indices`
        real(c_double), intent(in), dimension(n_values_elements), target :: values
            !! Input real array
        integer(c_int), intent(in), dimension(n_sorted_indices_elements), target :: sorted_indices
            !! Permutation index array
        integer(c_int), intent(in), target :: position
            !! Sorted position (1-based)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(out), target :: sorted_value
            !! 
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(values)
        M_CHECK_NON_NULL(n_values_elements)
        M_CHECK_NON_NULL(sorted_indices)
        M_CHECK_NON_NULL(n_sorted_indices_elements)
        M_CHECK_NON_NULL(position)
        M_CHECK_NON_NULL(sorted_value)
        sorted_value = get_sorted_value(&
            values = values,&
            sorted_indices = sorted_indices,&
            position = position,&
            ierr = ierr&
        )
    end subroutine get_sorted_value_c

    !> summary: C-wrapper for [[f42_binary_search_tree(module):build_bst_index(subroutine)]]
    !| Build the BST index by sorting indices using values in x.
    subroutine build_bst_index_c(&
            values,&
            num_values,&
            sorted_indices,&
            tmp_left_stack,&
            tmp_right_stack,&
            ierr&
            ) bind(C, name="build_bst_index_c")
        use f42_binary_search_tree, only: build_bst_index
        use f42_binary_search_tree
        integer(c_int), intent(in), target :: num_values
            !! Number of elements in values array
        real(c_double), intent(in), dimension(num_values), target :: values
            !! Input real array to be indexed
        integer(c_int), intent(out), dimension(num_values), target :: sorted_indices
            !! Output permutation index
        integer(c_int), intent(out), dimension(num_values), target :: tmp_left_stack
            !! Manual stack for left indices
        integer(c_int), intent(out), dimension(num_values), target :: tmp_right_stack
            !! Manual stack for right indices
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(values)
        M_CHECK_NON_NULL(num_values)
        M_CHECK_NON_NULL(sorted_indices)
        M_CHECK_NON_NULL(tmp_left_stack)
        M_CHECK_NON_NULL(tmp_right_stack)
        call build_bst_index(&
            values = values,&
            num_values = num_values,&
            sorted_indices = sorted_indices,&
            tmp_left_stack = tmp_left_stack,&
            tmp_right_stack = tmp_right_stack,&
            ierr = ierr&
        )
    end subroutine build_bst_index_c

    !> summary: C-wrapper for [[f42_binary_search_tree(module):bst_range_query(subroutine)]]
    !| Perform a 1D range query over the sorted index.
    subroutine bst_range_query_c(&
            values,&
            sorted_indices,&
            num_values,&
            lower_bound,&
            upper_bound,&
            output_indices,&
            num_matches,&
            ierr&
            ) bind(C, name="bst_range_query_c")
        use f42_binary_search_tree, only: bst_range_query
        use f42_binary_search_tree
        integer(c_int), intent(in), target :: num_values
            !! Number of elements
        real(c_double), intent(in), dimension(num_values), target :: values
            !! Input real array
        integer(c_int), intent(in), dimension(num_values), target :: sorted_indices
            !! Permutation index array (sorted)
        real(c_double), intent(in), target :: lower_bound
            !! Lower bound of range (inclusive)
        real(c_double), intent(in), target :: upper_bound
            !! Upper bound of range (inclusive)
        integer(c_int), intent(out), dimension(num_values), target :: output_indices
            !! Output array of matching indices
            !! The first `num_matches` elements will hold the results.
        integer(c_int), intent(out), target :: num_matches
            !! Number of matches found
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(values)
        M_CHECK_NON_NULL(sorted_indices)
        M_CHECK_NON_NULL(num_values)
        M_CHECK_NON_NULL(lower_bound)
        M_CHECK_NON_NULL(upper_bound)
        M_CHECK_NON_NULL(output_indices)
        M_CHECK_NON_NULL(num_matches)
        call bst_range_query(&
            values = values,&
            sorted_indices = sorted_indices,&
            num_values = num_values,&
            lower_bound = lower_bound,&
            upper_bound = upper_bound,&
            output_indices = output_indices,&
            num_matches = num_matches,&
            ierr = ierr&
        )
    end subroutine bst_range_query_c

end module f42_binary_search_tree_c
#endif