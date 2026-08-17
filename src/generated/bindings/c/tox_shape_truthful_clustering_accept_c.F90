#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_accept(module)]]
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
module tox_shape_truthful_clustering_accept_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: accept_ensemble_c
    public :: accept_ensemble_expert_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept(module):accept_ensemble(subroutine)]]
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
    subroutine accept_ensemble_c(&
            n_dimensions,&
            o,&
            U_first,&
            d_first,&
            U_history,&
            d_history,&
            history_len,&
            G_t,&
            normal_error_t,&
            U_tp1,&
            d_tp1,&
            G_tp1,&
            normal_error_tp1,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
            is_accepted,&
            ierr&
        ) bind(C, name="accept_ensemble_c")
        use tox_shape_truthful_clustering_accept, only: accept_ensemble

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth; sizes U_history/d_history below
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_first
            !! Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
            !! `misc/mod_STC.md`, "Output"
        integer(c_int), intent(in), target :: d_first
            !! Ensemble's intrinsic dimension at its bootstrap iteration
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions, n_dimensions, o), intent(in), target :: U_history
            !! Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
            !! valid (see `history_len`)
        integer(c_int), dimension(o), intent(in), target :: d_history
            !! Trailing intrinsic dimensions, one per column of U_history; only entries
            !! 1:history_len are valid
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: history_len
            !! Number of valid columns in U_history/d_history; column history_len is the most
            !! recently accepted iteration
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `o`.
        real(c_double), intent(in), target :: G_t
            !! Ensemble's spectral gap at the most recently accepted iteration
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: normal_error_t
            !! Ensemble's normal_error at the most recently accepted iteration, see
            !! `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
            !! noise at all in its normal directions -- unlike G_t, which is a ratio already
            !! protected by its own +epsilon denominator (see `observable`) and so is
            !! required to be strictly positive.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_tp1
            !! Candidate ensemble's tangent+normal basis
        integer(c_int), intent(in), target :: d_tp1
            !! Candidate ensemble's intrinsic dimension
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), target :: G_tp1
            !! Candidate ensemble's spectral gap
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: normal_error_tp1
            !! Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(out), target :: is_accepted
            !! .true. if all four acceptance criteria are satisfied
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(d_first)
        M_CHECK_NON_NULL(history_len)
        M_CHECK_NON_NULL(G_t)
        M_CHECK_NON_NULL(normal_error_t)
        M_CHECK_NON_NULL(d_tp1)
        M_CHECK_NON_NULL(G_tp1)
        M_CHECK_NON_NULL(normal_error_tp1)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(is_accepted)
        M_CHECK_ARRAY_NON_NULL(U_first, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(U_history, n_dimensions * n_dimensions * o)
        M_CHECK_ARRAY_NON_NULL(d_history, o)
        M_CHECK_ARRAY_NON_NULL(U_tp1, n_dimensions * n_dimensions)

        call accept_ensemble(&
            n_dimensions = n_dimensions,&
            o = o,&
            U_first = U_first,&
            d_first = d_first,&
            U_history = U_history,&
            d_history = d_history,&
            history_len = history_len,&
            G_t = G_t,&
            normal_error_t = normal_error_t,&
            U_tp1 = U_tp1,&
            d_tp1 = d_tp1,&
            G_tp1 = G_tp1,&
            normal_error_tp1 = normal_error_tp1,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
            is_accepted = is_accepted,&
            ierr = ierr&
        )
    end subroutine accept_ensemble_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept(module):accept_ensemble_expert(subroutine)]]
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
    subroutine accept_ensemble_expert_c(&
            n_dimensions,&
            o,&
            U_first,&
            d_first,&
            U_history,&
            d_history,&
            history_len,&
            G_t,&
            normal_error_t,&
            U_tp1,&
            d_tp1,&
            G_tp1,&
            normal_error_tp1,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
            lwork,&
            tmp_m,&
            tmp_s,&
            tmp_work,&
            is_accepted,&
            ierr&
        ) bind(C, name="accept_ensemble_expert_c")
        use tox_shape_truthful_clustering_accept, only: accept_ensemble_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth; sizes U_history/d_history below
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_accept_impl(module):tox_stc_accept_ensemble_svd_workspace]].
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_first
            !! Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
            !! `misc/mod_STC.md`, "Output"
        integer(c_int), intent(in), target :: d_first
            !! Ensemble's intrinsic dimension at its bootstrap iteration
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions, n_dimensions, o), intent(in), target :: U_history
            !! Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
            !! valid (see `history_len`)
        integer(c_int), dimension(o), intent(in), target :: d_history
            !! Trailing intrinsic dimensions, one per column of U_history; only entries
            !! 1:history_len are valid
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: history_len
            !! Number of valid columns in U_history/d_history; column history_len is the most
            !! recently accepted iteration
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `o`.
        real(c_double), intent(in), target :: G_t
            !! Ensemble's spectral gap at the most recently accepted iteration
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: normal_error_t
            !! Ensemble's normal_error at the most recently accepted iteration, see
            !! `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
            !! noise at all in its normal directions -- unlike G_t, which is a ratio already
            !! protected by its own +epsilon denominator (see `observable`) and so is
            !! required to be strictly positive.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_tp1
            !! Candidate ensemble's tangent+normal basis
        integer(c_int), intent(in), target :: d_tp1
            !! Candidate ensemble's intrinsic dimension
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), target :: G_tp1
            !! Candidate ensemble's spectral gap
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: normal_error_tp1
            !! Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(out), target :: tmp_m
            !! Workspace: M = U_r(:,1:d)^T U_tp1(:,1:d) for whichever reference is being compared
        real(c_double), dimension(n_dimensions), intent(out), target :: tmp_s
            !! Workspace: singular values of tmp_m (= cos of the principal angles)
        real(c_double), dimension(lwork), intent(out), target :: tmp_work
            !! Workspace: LAPACK dgesvd scratch
        logical(c_bool), intent(out), target :: is_accepted
            !! .true. if all four acceptance criteria are satisfied
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(d_first)
        M_CHECK_NON_NULL(history_len)
        M_CHECK_NON_NULL(G_t)
        M_CHECK_NON_NULL(normal_error_t)
        M_CHECK_NON_NULL(d_tp1)
        M_CHECK_NON_NULL(G_tp1)
        M_CHECK_NON_NULL(normal_error_tp1)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(lwork)
        M_CHECK_NON_NULL(is_accepted)
        M_CHECK_ARRAY_NON_NULL(U_first, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(U_history, n_dimensions * n_dimensions * o)
        M_CHECK_ARRAY_NON_NULL(d_history, o)
        M_CHECK_ARRAY_NON_NULL(U_tp1, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_m, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_s, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_work, lwork)

        call accept_ensemble_expert(&
            n_dimensions = n_dimensions,&
            o = o,&
            U_first = U_first,&
            d_first = d_first,&
            U_history = U_history,&
            d_history = d_history,&
            history_len = history_len,&
            G_t = G_t,&
            normal_error_t = normal_error_t,&
            U_tp1 = U_tp1,&
            d_tp1 = d_tp1,&
            G_tp1 = G_tp1,&
            normal_error_tp1 = normal_error_tp1,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
            lwork = lwork,&
            tmp_m = tmp_m,&
            tmp_s = tmp_s,&
            tmp_work = tmp_work,&
            is_accepted = is_accepted,&
            ierr = ierr&
        )
    end subroutine accept_ensemble_expert_c

end module tox_shape_truthful_clustering_accept_c
#endif
