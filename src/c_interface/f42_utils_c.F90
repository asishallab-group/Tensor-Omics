#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[f42_utils(module)]]
!| Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.
module f42_utils_c
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

    !> summary: C-wrapper for [[f42_utils(module):which(subroutine)]]
    !| Finds the indices of the true values in a logical mask.
    subroutine which_c(&
            mask,&
            n_mask_elements,&
            n,&
            idx_out,&
            n_idx_out_elements,&
            m_max,&
            m_out,&
            ierr&
            ) bind(C, name="which_c")
        use f42_utils, only: which
        use f42_utils
        integer(c_int), intent(in), target :: n_mask_elements
            !! Size of the 1. dimension/extent of `mask`
        integer(c_int), intent(in), target :: n_idx_out_elements
            !! Size of the 1. dimension/extent of `idx_out`
        integer(c_int), intent(in), dimension(n_mask_elements), target :: mask
            !! Logical array of size n.
        integer(c_int), intent(in), target :: n
            !! Size of the mask.
        integer(c_int), intent(out), dimension(n_idx_out_elements), target :: idx_out
            !! Integer array to store the indices of true values.
        integer(c_int), intent(in), target :: m_max
            !! Maximum size of idx_out.
        integer(c_int), intent(out), target :: m_out
            !! Actual size of idx_out (number of true values found).
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0=ok, 201=invalid input, 202=empty input
        logical, allocatable, dimension(:) :: mask_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(mask)
        M_CHECK_NON_NULL(n_mask_elements)
        M_CHECK_NON_NULL(n)
        M_CHECK_NON_NULL(idx_out)
        M_CHECK_NON_NULL(n_idx_out_elements)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(m_out)
        M_ALLOCATE(mask_f(n_mask_elements))
        call c_int_as_logical(mask, mask_f)
        call which(&
            mask = mask_f,&
            n = n,&
            idx_out = idx_out,&
            m_max = m_max,&
            m_out = m_out,&
            ierr = ierr&
        )
    end subroutine which_c

    !> summary: C-wrapper for [[f42_utils(module):loess_smooth_2d(subroutine)]]
    !| Performs LOESS smoothing on a set of data points.
    !| Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters.
    !| The user must pre-filter data and provide only valid indices in indices_used.
    subroutine loess_smooth_2d_c(&
            n_total,&
            n_target,&
            x_ref,&
            y_ref,&
            indices_used,&
            n_used,&
            x_query,&
            kernel_sigma,&
            kernel_cutoff,&
            y_out,&
            ierr&
            ) bind(C, name="loess_smooth_2d_c")
        use f42_utils, only: loess_smooth_2d
        use f42_utils
        integer(c_int), intent(in), target :: n_total
            !! Total number of reference points.
        integer(c_int), intent(in), target :: n_target
            !! Number of target points to smooth.
        integer(c_int), intent(in), target :: n_used
            !! Number of indices actually used for smoothing.
        real(c_double), intent(in), dimension(n_total), target :: x_ref
            !! Reference x-coordinates.
        real(c_double), intent(in), dimension(n_total), target :: y_ref
            !! Reference y-coordinates (length n_total).
        integer(c_int), intent(in), dimension(n_used), target :: indices_used
            !! Indices of reference points used for smoothing (only valid indices).
        real(c_double), intent(in), dimension(n_target), target :: x_query
            !! Target x-coordinates to smooth.
        real(c_double), intent(in), target :: kernel_sigma
            !! Bandwidth parameter for the kernel.
        real(c_double), intent(in), target :: kernel_cutoff
            !! Cutoff for the kernel.
        real(c_double), intent(out), dimension(n_target), target :: y_out
            !! Output smoothed values (length n_target).
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0=ok, 201=invalid input, 202=empty input
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_total)
        M_CHECK_NON_NULL(n_target)
        M_CHECK_NON_NULL(x_ref)
        M_CHECK_NON_NULL(y_ref)
        M_CHECK_NON_NULL(indices_used)
        M_CHECK_NON_NULL(n_used)
        M_CHECK_NON_NULL(x_query)
        M_CHECK_NON_NULL(kernel_sigma)
        M_CHECK_NON_NULL(kernel_cutoff)
        M_CHECK_NON_NULL(y_out)
        call loess_smooth_2d(&
            n_total = n_total,&
            n_target = n_target,&
            x_ref = x_ref,&
            y_ref = y_ref,&
            indices_used = indices_used,&
            n_used = n_used,&
            x_query = x_query,&
            kernel_sigma = kernel_sigma,&
            kernel_cutoff = kernel_cutoff,&
            y_out = y_out,&
            ierr = ierr&
        )
    end subroutine loess_smooth_2d_c

    !> summary: C-wrapper for [[f42_utils(module):compute_edf(subroutine)]]
    !| Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| Assumes perm is already sorted by values[perm]. Caller controls sorting algorithm.
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    subroutine compute_edf_expert_c(&
            values,&
            n_values,&
            perm,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
            ) bind(C, name="compute_edf_expert_c")
        use f42_utils, only: compute_edf
        use f42_utils
        integer(c_int), intent(in), target :: n_values
            !! Number of values in the input array.
        real(c_double), intent(in), dimension(n_values), target :: values
            !! Array of observed data values (e.g., contributions or spikes).
        integer(c_int), intent(in), dimension(n_values), target :: perm
            !! Pre-sorted permutation indices (must be sorted by values[perm]).
        real(c_double), intent(out), dimension(n_values), target :: unique_values
            !! Sorted unique data values.
        real(c_double), intent(out), dimension(n_values), target :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
        integer(c_int), intent(out), target :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0=ok, 201=invalid input, 202=empty input
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(values)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(perm)
        M_CHECK_NON_NULL(unique_values)
        M_CHECK_NON_NULL(cdf_values)
        M_CHECK_NON_NULL(n_unique)
        call compute_edf(&
            values = values,&
            n_values = n_values,&
            perm = perm,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique,&
            ierr = ierr&
        )
    end subroutine compute_edf_expert_c

    !> summary: C-wrapper for [[f42_utils(module):compute_edf_alloc(subroutine)]]
    !| Helper routine that sorts and calls compute_edf.
    !| Allocates workspace internally and performs sorting before computing EDF.
    !| Use this for convenience; use compute_edf directly for custom sorting.
    subroutine compute_edf_c(&
            values,&
            n_values,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
            ) bind(C, name="compute_edf_c")
        use f42_utils, only: compute_edf_alloc
        use f42_utils
        integer(c_int), intent(in), target :: n_values
            !! Number of values in the input array.
        real(c_double), intent(in), dimension(n_values), target :: values
            !! Array of observed data values (e.g., contributions or spikes).
        real(c_double), intent(out), dimension(n_values), target :: unique_values
            !! Sorted unique data values.
        real(c_double), intent(out), dimension(n_values), target :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
        integer(c_int), intent(out), target :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(c_int), intent(out), target :: ierr
            !! Error code: 0=ok, 201=invalid input, 202=empty input
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(values)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(unique_values)
        M_CHECK_NON_NULL(cdf_values)
        M_CHECK_NON_NULL(n_unique)
        call compute_edf_alloc(&
            values = values,&
            n_values = n_values,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique,&
            ierr = ierr&
        )
    end subroutine compute_edf_c

    !> summary: C-wrapper for [[f42_utils(module):compute_empirical_p_values(subroutine)]]
    !| Calculate empirical p-values for scaled expression distances (RDI).
    !| 
    !| Implements:
    !| P(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !| 
    !| Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
    !| 
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    subroutine compute_empirical_p_values_c(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            perm,&
            p_values,&
            c_const,&
            ierr&
            ) bind(C, name="compute_empirical_p_values_c")
        use f42_utils, only: compute_empirical_p_values
        use f42_utils
        integer(c_int), intent(in), target :: n_genes
            !! 
        real(c_double), intent(in), dimension(n_genes), target :: rdi
            !! Number of genes being processed.
        real(c_double), intent(in), dimension(n_genes), target :: sorted_rdi
            !! empirical distribution D
        integer(c_int), intent(in), dimension(n_genes), target :: perm
            !! Constant used in the computation, typically 1
        real(c_double), intent(out), dimension(n_genes), target :: p_values
            !! empirical distribution D with non negative values
        real(c_double), intent(in), target :: c_const
            !! Output array to store the computed p-values for each gene.
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(rdi)
        M_CHECK_NON_NULL(sorted_rdi)
        M_CHECK_NON_NULL(perm)
        M_CHECK_NON_NULL(p_values)
        M_CHECK_NON_NULL(c_const)
        call compute_empirical_p_values(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            p_values = p_values,&
            c_const = c_const&
        )
    end subroutine compute_empirical_p_values_c

end module f42_utils_c
#endif