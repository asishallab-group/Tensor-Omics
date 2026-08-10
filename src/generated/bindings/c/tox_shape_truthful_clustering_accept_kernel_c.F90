#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_accept_kernel(module)]]
!| # Shape Truthful Clustering (STC): Accept
!|
!| `accept_ensemble`: whether a grown ensemble at t+1 is still compatible with its own growth
!| trajectory, judged by four criteria -- tangent-space drift (chordal distance, compared
!| against a reference set: the bootstrap iteration plus the trailing o-window, not just the
!| immediately preceding iteration), change in intrinsic dimension (against both the bootstrap
!| iteration and the immediately preceding one), relative change in spectral gap, and relative
!| change in residual (RMSE), both against the immediately preceding iteration only. See
!| `misc/mod_STC.md`, SKG `accept_ensemble`, for the full algorithm definition and the
!| "no cumulative-rotation budget" rationale for comparing against a reference set rather than
!| a single previous state. This compares the SAME ensemble across one growth step -- not two
!| different ensembles/anchors at a possible junction -- so, unlike
!| `misc/STC_for_LoManLe.md` section 4's explicit "angle never gates a junction" rule, a
!| tangent-space-drift mismatch here legitimately contributes to rejection.
module tox_shape_truthful_clustering_accept_kernel_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: tox_stc_accept_ensemble_svd_workspace_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept_kernel(module):tox_stc_accept_ensemble_svd_workspace(subroutine)]]
    !| Sized for the worst case across every reference comparison accept_ensemble performs (up
    !| to n_dimensions-square, since a comparison's actual shared rank is always
    !| <= n_dimensions): the documented minimum-workspace formula for a square M=N=n_dimensions
    !| input with JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*n_dimensions). A
    !| larger-than-required LWORK is always safe per LAPACK's own convention, so one size,
    !| computed once, serves every one of accept_ensemble's (up to o+1) comparisons.
    subroutine tox_stc_accept_ensemble_svd_workspace_c(&
            n_dimensions,&
            lwork,&
            ierr&
        ) bind(C, name="tox_stc_accept_ensemble_svd_workspace_c")
        use tox_shape_truthful_clustering_accept_kernel, only: tox_stc_accept_ensemble_svd_workspace

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(out), target :: lwork
            !! Recommended size of the real LAPACK workspace
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(lwork)

        call tox_stc_accept_ensemble_svd_workspace(&
            n_dimensions = n_dimensions,&
            lwork = lwork&
        )
    end subroutine tox_stc_accept_ensemble_svd_workspace_c

end module tox_shape_truthful_clustering_accept_kernel_c
#endif
