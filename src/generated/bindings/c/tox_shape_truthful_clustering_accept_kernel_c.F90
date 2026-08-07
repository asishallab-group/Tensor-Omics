#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_accept_kernel(module)]]
!| # Shape Truthful Clustering (STC): Accept
!|
!| `accept_ensemble`: whether a grown ensemble at t+1 is still compatible with its own state
!| at t, judged by principal angle between tangent bases, change in intrinsic dimension, and
!| relative change in spectral gap. See `misc/mod_STC.md`, SKG `accept_ensemble`, for the
!| full algorithm definition. This compares the SAME ensemble across one growth step -- not
!| two different ensembles/anchors at a possible junction -- so, unlike
!| `misc/STC_for_LoManLe.md` section 4's explicit "angle never gates a junction" rule, a
!| principal-angle mismatch here legitimately contributes to rejection.
module tox_shape_truthful_clustering_accept_kernel_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: tox_stc_accept_ensemble_svd_workspace_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept_kernel(module):tox_stc_accept_ensemble_svd_workspace(subroutine)]]
    !| The documented minimum-workspace formula for a square M=N=min(d_t,d_tp1) input with
    !| JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*min(M,N)).
    subroutine tox_stc_accept_ensemble_svd_workspace_c(&
            d_t,&
            d_tp1,&
            lwork,&
            ierr&
        ) bind(C, name="tox_stc_accept_ensemble_svd_workspace_c")
        use tox_shape_truthful_clustering_accept_kernel, only: tox_stc_accept_ensemble_svd_workspace

        integer(c_int), intent(in), target :: d_t
            !! Ensemble's intrinsic dimension at t
        integer(c_int), intent(in), target :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
        integer(c_int), intent(out), target :: lwork
            !! Recommended size of the real LAPACK workspace
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(d_t)
        M_CHECK_NON_NULL(d_tp1)
        M_CHECK_NON_NULL(lwork)

        call tox_stc_accept_ensemble_svd_workspace(&
            d_t = d_t,&
            d_tp1 = d_tp1,&
            lwork = lwork&
        )
    end subroutine tox_stc_accept_ensemble_svd_workspace_c

end module tox_shape_truthful_clustering_accept_kernel_c
#endif
