#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_preprocessing_kernel(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| Kernels for the pipeline that turns expression vectors into neighborhood residuals, the
!| input of the JCT based data integration. The generator turns the `*_kernel` procedures
!| into the validating and allocating wrappers in module tox_data_integration_preprocessing;
!| `calc_neighborhood_size` is a recommend routine and is exported as it stands.
module tox_data_integration_preprocessing_kernel_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_neighborhood_size_c

contains

    !> summary: C-wrapper for [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]]
    !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower
    !| due to few genes with non-NaN mean.
    subroutine calc_neighborhood_size_c(&
            n_pool,&
            n_points,&
            n_genes_S,&
            mean_S,&
            desired_size,&
            n_neighbors,&
            ierr&
        ) bind(C, name="calc_neighborhood_size_c")
        use tox_data_integration_preprocessing_kernel, only: calc_neighborhood_size

        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_pool
            !! Total number of pooled mean-expression values across both studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: desired_size
            !! Optional desired neighborhood size
            !! The default value is `1000`.
        integer(c_int), intent(out), target :: n_neighbors
            !! Calculated neighborhood size
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(desired_size)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)

        n_neighbors = calc_neighborhood_size(&
            n_pool = n_pool,&
            n_points = n_points,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            desired_size = desired_size&
        )
    end subroutine calc_neighborhood_size_c

end module tox_data_integration_preprocessing_kernel_c
#endif
