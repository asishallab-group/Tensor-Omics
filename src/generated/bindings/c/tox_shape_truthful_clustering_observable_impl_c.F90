#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_observable_impl(module)]]
!| # Shape Truthful Clustering (STC): Observable
!|
!| `observable`: the tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble,
!| obtained from the economy-mode singular value decomposition (LAPACK `dgesdd`) of its
!| centered member vectors -- never an eigendecomposition of an explicitly formed
!| covariance matrix (see `misc/mod_STC.md`, "Numerical Linear Algebra"). `normal_error` and
!| `tangent_scales` are simple, dependency-free reductions over the eigenvalues `observable`
!| computes. See `misc/mod_STC.md`, SKG `observable`/`normal_error`/`tangent_scales`, for the
!| full algorithm definitions.
module tox_shape_truthful_clustering_observable_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: tox_stc_observable_svd_workspace_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable_impl(module):tox_stc_observable_svd_workspace(subroutine)]]
    !| The documented minimum-workspace formula for JOBZ='S' (see `man dgesdd`):
    !| LWORK >= 4*min(M,N)**2 + 7*min(M,N), IWORK size = 8*min(M,N), where M=n_dimensions and
    !| N=n_selected_member.
    subroutine tox_stc_observable_svd_workspace_c(&
            n_dimensions,&
            n_selected_member,&
            lwork,&
            iwork_size,&
            ierr&
        ) bind(C, name="tox_stc_observable_svd_workspace_c")
        use tox_shape_truthful_clustering_observable_impl, only: tox_stc_observable_svd_workspace

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_selected_member
            !! Number of selected ensemble members
        integer(c_int), intent(out), target :: lwork
            !! Recommended size of the real LAPACK workspace
        integer(c_int), intent(out), target :: iwork_size
            !! Recommended size of the integer LAPACK workspace
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_selected_member)
        M_CHECK_NON_NULL(lwork)
        M_CHECK_NON_NULL(iwork_size)

        call tox_stc_observable_svd_workspace(&
            n_dimensions = n_dimensions,&
            n_selected_member = n_selected_member,&
            lwork = lwork,&
            iwork_size = iwork_size&
        )
    end subroutine tox_stc_observable_svd_workspace_c

end module tox_shape_truthful_clustering_observable_impl_c
#endif
