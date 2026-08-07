#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Accept
!|
!| `accept_ensemble`: whether a grown ensemble at t+1 is still compatible with its own state
!| at t, judged by principal angle between tangent bases, change in intrinsic dimension, and
!| relative change in spectral gap. See `misc/mod_STC.md`, SKG `accept_ensemble`, for the
!| full algorithm definition. This compares the SAME ensemble across one growth step -- not
!| two different ensembles/anchors at a possible junction -- so, unlike
!| `misc/STC_for_LoManLe.md` section 4's explicit "angle never gates a junction" rule, a
!| principal-angle mismatch here legitimately contributes to rejection.
module tox_shape_truthful_clustering_accept_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err_once, ERR_INTERNAL
    use f42_utils, only: above
    M_IMPLICIT_NONE

    interface
        ! Declared pure for the same reason as dgesdd in
        ! tox_shape_truthful_clustering_observable_kernel. `u`/`vt` are assumed-size (`*`)
        ! since this module's only call site uses JOBU='N', JOBVT='N', for which LAPACK
        ! documents them as "not referenced" -- their exact shape is irrelevant.
        pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import :: int32, real64
            character,      intent(in)    :: jobu, jobvt
            integer(int32), intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(real64),   intent(inout) :: a(lda, n)
            real(real64),   intent(out)   :: s(min(m, n))
            real(real64),   intent(out)   :: u(ldu, *)
            real(real64),   intent(out)   :: vt(ldvt, *)
            real(real64),   intent(out)   :: work(lwork)
            integer(int32), intent(out)   :: info
        end subroutine dgesvd
    end interface

    private
    public :: accept_ensemble_kernel
    public :: tox_stc_accept_ensemble_svd_workspace

contains

    !> M_EXPORT_C
    !| summary: Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVD
    !| AUTHOR_ASIS_HALLAB
    !| The documented minimum-workspace formula for a square M=N=min(d_t,d_tp1) input with
    !| JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*min(M,N)).
    pure subroutine tox_stc_accept_ensemble_svd_workspace(d_t, d_tp1, lwork)
        integer(int32), intent(in) :: d_t
            !! Ensemble's intrinsic dimension at t
        integer(int32), intent(in) :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
        integer(int32), intent(out) :: lwork
            !! Recommended size of the real LAPACK workspace

        lwork = max(1_int32, 5_int32*min(d_t, d_tp1))

    end subroutine tox_stc_accept_ensemble_svd_workspace

    !> summary: Whether a grown ensemble at t+1 is still compatible with its own state at t
    !| AUTHOR_ASIS_HALLAB
    !| Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
    !| bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
    !| cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
    !| dimension itself changed, the two tangent bases don't share a common dimension to
    !| compare angles over at all, so this criterion is skipped (no SVD is computed) and
    !| treated as vacuously satisfied; criterion (2) is what actually judges whether that
    !| change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
    !| `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
    !| check could foresee.
    pure subroutine accept_ensemble_kernel(n_dimensions, U_t, d_t, G_t, U_tp1, d_tp1, G_tp1, &
                                           alpha_max, d_max, G_max, lwork, &
                                           tmp_m, tmp_s, tmp_work, &
                                           is_accepted, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        real(real64), intent(in) :: U_t(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t, see `observable`
        integer(int32), intent(in) :: d_t
            !! Ensemble's intrinsic dimension at t
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at t
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: U_tp1(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at t+1
        integer(int32), intent(in) :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: G_tp1
            !! Ensemble's spectral gap at t+1
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(2.0_real64 * atan(1.0_real64))
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! DM_OUTPUT_FROM(lwork, tox_stc_accept_ensemble_svd_workspace, tox_shape_truthful_clustering_accept_kernel, AUTO)
        real(real64), intent(out) :: tmp_m(min(d_t, d_tp1), min(d_t, d_tp1))
            !! Workspace: M = U_t(:,1:d)^T U_tp1(:,1:d)
        real(real64), intent(out) :: tmp_s(min(d_t, d_tp1))
            !! Workspace: singular values of M (= cos of the principal angles)
        real(real64), intent(out) :: tmp_work(lwork)
            !! Workspace: LAPACK dgesvd scratch
        logical, intent(out) :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        real(real64)   :: u_dummy(1, 1), vt_dummy(1, 1), cos_alpha, alpha_i
        logical        :: angle_ok, d_diff_ok, g_ratio_ok
        integer(int32) :: d_common, i, info

        call set_ok(ierr)

        d_diff_ok  = abs(d_tp1 - d_t) <= d_max
        g_ratio_ok = abs(log(G_tp1/G_t)) <= G_max

        angle_ok = .true.
        if (d_t == d_tp1) then
            d_common = d_t
            if (d_common > 0) then
                tmp_m = matmul(transpose(U_t(:, 1:d_common)), U_tp1(:, 1:d_common))
                call dgesvd('N', 'N', d_common, d_common, tmp_m, d_common, tmp_s, &
                           u_dummy, 1, vt_dummy, 1, tmp_work, lwork, info)
                if (info /= 0) then
                    call set_err_once(ierr, ERR_INTERNAL)
                    return
                end if
                do i = 1, d_common
                    cos_alpha = max(-1.0_real64, min(1.0_real64, tmp_s(i)))
                    alpha_i   = acos(cos_alpha)
                    if (alpha_i > alpha_max) angle_ok = .false.
                end do
            end if
        end if

        is_accepted = angle_ok .and. d_diff_ok .and. g_ratio_ok

    end subroutine accept_ensemble_kernel

end module tox_shape_truthful_clustering_accept_kernel
