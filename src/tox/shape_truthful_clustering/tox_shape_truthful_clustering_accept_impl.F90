#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Accept
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
module tox_shape_truthful_clustering_accept_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, ERR_INTERNAL
    use f42_math_impl, only: above
    M_IMPLICIT_NONE

    interface
        ! Declared pure for the same reason as dgesdd in
        ! tox_shape_truthful_clustering_observable_impl. `u`/`vt` are assumed-size (`*`)
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
    public :: accept_ensemble_impl
    public :: tox_stc_accept_ensemble_svd_workspace
    public :: stc_chordal_distance

contains

    !> M_EXPORT_C
    !| summary: Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVDs
    !| AUTHOR_ASIS_HALLAB
    !| Sized for the worst case across every reference comparison accept_ensemble performs (up
    !| to n_dimensions-square, since a comparison's actual shared rank is always
    !| <= n_dimensions): the documented minimum-workspace formula for a square M=N=n_dimensions
    !| input with JOBU='N', JOBVT='N' (see `man dgesvd`): LWORK >= max(1, 5*n_dimensions). A
    !| larger-than-required LWORK is always safe per LAPACK's own convention, so one size,
    !| computed once, serves every one of accept_ensemble's (up to o+1) comparisons.
    pure subroutine tox_stc_accept_ensemble_svd_workspace(n_dimensions, lwork)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(out) :: lwork
            !! Recommended size of the real LAPACK workspace

        lwork = max(1_int32, 5_int32*n_dimensions)

    end subroutine tox_stc_accept_ensemble_svd_workspace

    !> summary: Whether a grown ensemble at t+1 is still compatible with its own growth trajectory
    !| AUTHOR_ASIS_HALLAB
    !| Four criteria, all must hold: (1) tangent-space drift -- for every reference basis in
    !| {U_first} union {U_history(:,:,1:history_len)} that shares d_tp1's rank, the chordal
    !| distance to U_tp1 (via `dgesvd` on M = U_r(:,1:d)^T U_tp1(:,1:d), whose singular values
    !| are cos(theta_i) directly) must not exceed
    !| chordal_dist_max_as_prcnt_of_range * sqrt(d_tp1); references with a different d are
    !| skipped for this criterion (no shared dimension to compare angles over), judged instead
    !| by criterion (2). Only the candidate is compared against each reference, not every pair
    !| within the reference set itself -- sufficient, not a shortcut, since every reference
    !| already passed this same criterion against every other reference present in the set at
    !| the growth step it was itself accepted (an inductive invariant FIFO eviction cannot
    !| break), so this costs O(history_len) small SVDs per growth step, not O(history_len^2).
    !| (2) max(|d_tp1 - d_first|, |d_tp1 - d_history(history_len)|) <= d_max. (3)
    !| |log(G_tp1/G_t)| <= G_max, where G_t = the most recently accepted iteration's spectral
    !| gap. (4) |log(RMSE_tp1/RMSE_t)| <= RMSE_change_max, where RMSE = sqrt(normal_error) --
    !| see `normal_error` in the sibling observable kernel module; already free, no new SVD or
    !| storage. `ierr` is set only if a LAPACK SVD fails to converge -- not a condition any
    !| input check could foresee.
    pure subroutine accept_ensemble_impl(n_dimensions, o, &
                                           U_first, d_first, &
                                           U_history, d_history, history_len, &
                                           G_t, normal_error_t, &
                                           U_tp1, d_tp1, G_tp1, normal_error_tp1, &
                                           chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, &
                                           lwork, tmp_m, tmp_s, tmp_work, &
                                           is_accepted, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth; sizes U_history/d_history below
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: U_first(n_dimensions, n_dimensions)
            !! Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
            !! `misc/mod_STC.md`, "Output"
        integer(int32), intent(in) :: d_first
            !! Ensemble's intrinsic dimension at its bootstrap iteration
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: U_history(n_dimensions, n_dimensions, o)
            !! Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
            !! valid (see `history_len`)
        integer(int32), intent(in) :: d_history(o)
            !! Trailing intrinsic dimensions, one per column of U_history; only entries
            !! 1:history_len are valid
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: history_len
            !! Number of valid columns in U_history/d_history; column history_len is the most
            !! recently accepted iteration
            !! DM_MIN(1_int32)
            !! DM_MAX(o)
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at the most recently accepted iteration
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: normal_error_t
            !! Ensemble's normal_error at the most recently accepted iteration, see
            !! `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
            !! noise at all in its normal directions -- unlike G_t, which is a ratio already
            !! protected by its own +epsilon denominator (see `observable`) and so is
            !! required to be strictly positive.
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: U_tp1(n_dimensions, n_dimensions)
            !! Candidate ensemble's tangent+normal basis
        integer(int32), intent(in) :: d_tp1
            !! Candidate ensemble's intrinsic dimension
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: G_tp1
            !! Candidate ensemble's spectral gap
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: normal_error_tp1
            !! Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|
            !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! DM_OUTPUT_FROM(lwork, tox_stc_accept_ensemble_svd_workspace, tox_shape_truthful_clustering_accept_impl, AUTO)
        real(real64), intent(out) :: tmp_m(n_dimensions, n_dimensions)
            !! Workspace: M = U_r(:,1:d)^T U_tp1(:,1:d) for whichever reference is being compared
        real(real64), intent(out) :: tmp_s(n_dimensions)
            !! Workspace: singular values of tmp_m (= cos of the principal angles)
        real(real64), intent(out) :: tmp_work(lwork)
            !! Workspace: LAPACK dgesvd scratch
        logical(c_bool), intent(out) :: is_accepted
            !! .true. if all four acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        real(real64)   :: rmse_t, rmse_tp1, chordal_distance, chordal_max
        logical(c_bool)        :: applicable, chordal_ok, d_diff_ok, g_ratio_ok, rmse_ratio_ok
        integer(int32) :: i

        call set_ok(ierr)

        chordal_max = chordal_dist_max_as_prcnt_of_range * sqrt(real(d_tp1, real64))
        chordal_ok  = .true.

        call stc_chordal_distance(n_dimensions, U_first, d_first, U_tp1, d_tp1, lwork, tmp_m, tmp_s, tmp_work, &
                                  applicable, chordal_distance, ierr)
        if (ierr /= 0) return
        if (applicable .and. chordal_distance > chordal_max) chordal_ok = .false.

        do i = 1, history_len
            call stc_chordal_distance(n_dimensions, U_history(:, :, i), d_history(i), U_tp1, d_tp1, &
                                      lwork, tmp_m, tmp_s, tmp_work, applicable, chordal_distance, ierr)
            if (ierr /= 0) return
            if (applicable .and. chordal_distance > chordal_max) chordal_ok = .false.
        end do

        d_diff_ok = max(abs(d_tp1 - d_first), abs(d_tp1 - d_history(history_len))) <= d_max

        g_ratio_ok = abs(log(G_tp1/G_t)) <= G_max

        ! +epsilon: normal_error, unlike G, is a raw sum of eigenvalues with no ratio-internal
        ! protection of its own, and can be genuinely, exactly zero for a perfectly
        ! flat/collinear ensemble -- without this, log(RMSE_tp1/RMSE_t) is log(0/0), undefined,
        ! the first time growth is ever perfectly noise-free. See misc/mod_STC.md, criterion (4).
        rmse_t        = sqrt(normal_error_t + epsilon(1.0_real64))
        rmse_tp1      = sqrt(normal_error_tp1 + epsilon(1.0_real64))
        rmse_ratio_ok = abs(log(rmse_tp1/rmse_t)) <= RMSE_change_max

        is_accepted = chordal_ok .and. d_diff_ok .and. g_ratio_ok .and. rmse_ratio_ok

    end subroutine accept_ensemble_impl

    !> Chordal distance (Edelman, Arias & Smith 1998; related to the Davis-Kahan sinTheta
    !| theorem) between two tangent bases of matching rank, or "not applicable" when their
    !| ranks differ -- no shared dimension to compare angles over, see `accept_ensemble`'s
    !| criterion (1). Not itself a kernel (no generator wrapper -- callers pass their own
    !| workspace directly): shared scratch-reuse plumbing for `accept_ensemble_impl`'s
    !| reference-set loop, exported `public` so `tox_stc_json`'s report-layer drift
    !| computations (see `misc/mod_STC.md`, "Ensemble Observable Plots") can reuse the exact
    !| same formula instead of re-deriving it. `ierr` is set only on a genuine LAPACK SVD
    !| non-convergence.
    pure subroutine stc_chordal_distance(n_dimensions, U_a, d_a, U_b, d_b, lwork, tmp_m, tmp_s, tmp_work, &
                                         applicable, chordal_distance, ierr)
        integer(int32), intent(in)    :: n_dimensions, d_a, d_b, lwork
        real(real64),   intent(in)    :: U_a(n_dimensions, n_dimensions)
        real(real64),   intent(in)    :: U_b(n_dimensions, n_dimensions)
        real(real64),   intent(out)   :: tmp_m(n_dimensions, n_dimensions)
        real(real64),   intent(out)   :: tmp_s(n_dimensions)
        real(real64),   intent(out)   :: tmp_work(lwork)
        logical(c_bool),        intent(out)   :: applicable
        real(real64),   intent(out)   :: chordal_distance
        integer(int32), intent(out)   :: ierr

        real(real64)   :: u_dummy(1, 1), vt_dummy(1, 1), cos_theta
        integer(int32) :: d_common, i, info

        call set_ok(ierr)
        chordal_distance = 0.0_real64
        applicable        = (d_a == d_b .and. d_a > 0)
        if (.not. applicable) return

        d_common = d_a
        tmp_m(1:d_common, 1:d_common) = matmul(transpose(U_a(:, 1:d_common)), U_b(:, 1:d_common))
        call dgesvd('N', 'N', d_common, d_common, tmp_m(1:d_common, 1:d_common), d_common, tmp_s(1:d_common), &
                   u_dummy, 1, vt_dummy, 1, tmp_work, lwork, info)
        if (info /= 0) then
            call set_err_once(ierr, ERR_INTERNAL)
            return
        end if

        do i = 1, d_common
            cos_theta         = max(-1.0_real64, min(1.0_real64, tmp_s(i)))
            chordal_distance = chordal_distance + (1.0_real64 - cos_theta**2)
        end do
        chordal_distance = sqrt(chordal_distance)

    end subroutine stc_chordal_distance

end module tox_shape_truthful_clustering_accept_impl
