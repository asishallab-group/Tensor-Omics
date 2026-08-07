#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_stats_impl(module)]]
!| Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_stats_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: loess_smooth_2d_c
    public :: compute_scaled_distance_quantile_c

contains

    !> summary: C-wrapper for [[f42_stats_impl(module):loess_smooth_2d(subroutine)]]
    !| Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
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
        use f42_stats_impl, only: loess_smooth_2d

        integer(c_int), intent(in), target :: n_total
            !! Total number of reference points.
        integer(c_int), intent(in), target :: n_target
            !! Number of target points to smooth.
        integer(c_int), intent(in), target :: n_used
            !! Number of indices actually used for smoothing.
        real(c_double), dimension(n_total), intent(in), target :: x_ref
            !! Reference x-coordinates.
        real(c_double), dimension(n_total), intent(in), target :: y_ref
            !! Reference y-coordinates (length n_total).
        integer(c_int), dimension(n_used), intent(in), target :: indices_used
            !! Indices of reference points used for smoothing (only valid indices).
        real(c_double), dimension(n_target), intent(in), target :: x_query
            !! Target x-coordinates to smooth.
        real(c_double), intent(in), target :: kernel_sigma
            !! Bandwidth parameter for the kernel.
        real(c_double), intent(in), target :: kernel_cutoff
            !! Cutoff for the kernel, not used if zero
        real(c_double), dimension(n_target), intent(out), target :: y_out
            !! Output smoothed values (length n_target).
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_total)
        M_CHECK_NON_NULL(n_target)
        M_CHECK_NON_NULL(n_used)
        M_CHECK_NON_NULL(kernel_sigma)
        M_CHECK_NON_NULL(kernel_cutoff)
        M_CHECK_ARRAY_NON_NULL(x_ref, n_total)
        M_CHECK_ARRAY_NON_NULL(y_ref, n_total)
        M_CHECK_ARRAY_NON_NULL(indices_used, n_used)
        M_CHECK_ARRAY_NON_NULL(x_query, n_target)
        M_CHECK_ARRAY_NON_NULL(y_out, n_target)

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

    !> summary: C-wrapper for [[f42_stats_impl(module):compute_scaled_distance_quantile(subroutine)]]
    !| This is NOT a null-hypothesis-testing p-value: each distance is compared against the
    !| observed distribution it was drawn from, not an independently generated null distribution.
    !| It instead measures how extreme an observed distance is relative to all observed distances.
    !|
    !| Implements:
    !| Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail quantile is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    subroutine compute_scaled_distance_quantile_c(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            perm,&
            quantile,&
            c_const,&
            ierr&
        ) bind(C, name="compute_scaled_distance_quantile_c")
        use f42_stats_impl, only: compute_scaled_distance_quantile

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes being processed.
        real(c_double), dimension(n_genes), intent(in), target :: rdi
            !! empirical distribution D
        real(c_double), dimension(n_genes), intent(in), target :: sorted_rdi
            !! empirical distribution D with non negative values
        integer(c_int), dimension(n_genes), intent(in), target :: perm
            !! Permutation array with sorted indices for sorted_rdi
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Output array to store the computed quantile for each gene.
        real(c_double), intent(in), target :: c_const
            !! Constant used in the computation, typically 1
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(c_const)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call compute_scaled_distance_quantile(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            perm = perm,&
            quantile = quantile,&
            c_const = c_const&
        )
    end subroutine compute_scaled_distance_quantile_c

end module f42_stats_impl_c
#endif
