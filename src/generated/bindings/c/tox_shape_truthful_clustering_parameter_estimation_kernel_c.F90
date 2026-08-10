#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_parameter_estimation_kernel(module)]]
!| # Shape Truthful Clustering (STC): Parameter Estimation
!|
!| A separate, optional pipeline step estimating near-optimal starting values for the crucial
!| parameters (`k_min`, `k_density`, `density_quantile`,
!| `chordal_dist_max_as_prcnt_of_range`, `G_max`, `d_max`) directly from the input data, at a
!| fraction of the cost of a grid search or a
!| resampling-based scheme: grow a handful of "estimator anchors" (EAs) into small local
!| neighborhoods using the same primitives the real pipeline already has
!| (`density_labels`, `observable`), then read the parameters off simple summary statistics of
!| that pass. See `misc/mod_STC.md`, "Estimate parameters from data", for the full algorithm
!| definition and the reasoning behind every design choice below.
module tox_shape_truthful_clustering_parameter_estimation_kernel_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: tox_stc_estimate_parameters_svd_workspace_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace(subroutine)]]
    !| Worst-case sizing for both `observable`'s dgesdd (an EA cloud can be as large as
    !| n_vectors) and the pairwise principal-angle dgesvd (shared rank at most n_dimensions),
    !| see `tox_stc_observable_svd_workspace` and `tox_stc_accept_ensemble_svd_workspace` for
    !| the individual formulas this combines.
    subroutine tox_stc_estimate_parameters_svd_workspace_c(&
            n_dimensions,&
            n_vectors,&
            lwork_observable,&
            iwork_size,&
            lwork_angle,&
            ierr&
        ) bind(C, name="tox_stc_estimate_parameters_svd_workspace_c")
        use tox_shape_truthful_clustering_parameter_estimation_kernel, only: tox_stc_estimate_parameters_svd_workspace

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(out), target :: lwork_observable
            !! Recommended size of observable's real LAPACK workspace (worst case)
        integer(c_int), intent(out), target :: iwork_size
            !! Recommended size of observable's integer LAPACK workspace (worst case)
        integer(c_int), intent(out), target :: lwork_angle
            !! Recommended size of the pairwise principal-angle LAPACK workspace (worst case)
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(lwork_observable)
        M_CHECK_NON_NULL(iwork_size)
        M_CHECK_NON_NULL(lwork_angle)

        call tox_stc_estimate_parameters_svd_workspace(&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            lwork_observable = lwork_observable,&
            iwork_size = iwork_size,&
            lwork_angle = lwork_angle&
        )
    end subroutine tox_stc_estimate_parameters_svd_workspace_c

end module tox_shape_truthful_clustering_parameter_estimation_kernel_c
#endif
