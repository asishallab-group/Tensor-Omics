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
!|
!| Generated from [[tox_shape_truthful_clustering_accept_impl(module)]]; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_accept
    use f42_safeguard
    use tox_shape_truthful_clustering_accept_impl, only: accept_ensemble_impl, tox_stc_accept_ensemble_svd_workspace
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: accept_ensemble
    public :: accept_ensemble_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_shape_truthful_clustering_accept_impl(module):accept_ensemble_impl]] needs, then calls it. The entry point to reach for first; see [[tox_shape_truthful_clustering_accept(module):accept_ensemble_expert]] to prepare it yourself.
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
    pure subroutine accept_ensemble(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth; sizes U_history/d_history below
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_first
            !! Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
            !! `misc/mod_STC.md`, "Output"
        integer(int32), intent(in) :: d_first
            !! Ensemble's intrinsic dimension at its bootstrap iteration
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions, n_dimensions, o), intent(in) :: U_history
            !! Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
            !! valid (see `history_len`)
        integer(int32), dimension(o), intent(in) :: d_history
            !! Trailing intrinsic dimensions, one per column of U_history; only entries
            !! 1:history_len are valid
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: history_len
            !! Number of valid columns in U_history/d_history; column history_len is the most
            !! recently accepted iteration
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `o`.
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at the most recently accepted iteration
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: normal_error_t
            !! Ensemble's normal_error at the most recently accepted iteration, see
            !! `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
            !! noise at all in its normal directions -- unlike G_t, which is a ratio already
            !! protected by its own +epsilon denominator (see `observable`) and so is
            !! required to be strictly positive.
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_tp1
            !! Candidate ensemble's tangent+normal basis
        integer(int32), intent(in) :: d_tp1
            !! Candidate ensemble's intrinsic dimension
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), intent(in) :: G_tp1
            !! Candidate ensemble's spectral gap
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: normal_error_tp1
            !! Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in) :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in) :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(out) :: is_accepted
            !! .true. if all four acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success
        integer(int32) :: lwork
        real(real64), dimension(:, :), allocatable :: tmp_m
        real(real64), dimension(:), allocatable :: tmp_s
        real(real64), dimension(:), allocatable :: tmp_work

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=1_int32)
        call validate_in_range_int(o, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(d_first, ierr, arg_pos=4_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_int(history_len, ierr, arg_pos=7_int32, min=1_int32, max=o)
        call validate_in_range_real(G_t, ierr, arg_pos=8_int32, min=above(0.0_real64))
        call validate_in_range_real(normal_error_t, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_in_range_int(d_tp1, ierr, arg_pos=11_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_tp1, ierr, arg_pos=12_int32, min=above(0.0_real64))
        call validate_in_range_real(normal_error_tp1, ierr, arg_pos=13_int32, min=0.0_real64)
        call validate_in_range_real(chordal_dist_max_as_prcnt_of_range, ierr, arg_pos=14_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(d_max, ierr, arg_pos=15_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=16_int32, min=0.0_real64)
        call validate_in_range_real(RMSE_change_max, ierr, arg_pos=17_int32, min=0.0_real64)
        call validate_all_in_range_real(U_first, n_dimensions * n_dimensions, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(U_history, n_dimensions * n_dimensions * o, ierr, arg_pos=5_int32)
        call validate_all_in_range_int(d_history, o, ierr, arg_pos=6_int32, min=0_int32, max=n_dimensions)
        call validate_all_in_range_real(U_tp1, n_dimensions * n_dimensions, ierr, arg_pos=10_int32)
        if (is_err(ierr)) return
#endif

        call tox_stc_accept_ensemble_svd_workspace(&
            n_dimensions = n_dimensions,&
            lwork = lwork&
        )
        M_ALLOCATE(tmp_m(n_dimensions, n_dimensions))
        M_ALLOCATE(tmp_s(n_dimensions))
        M_ALLOCATE(tmp_work(lwork))

        call accept_ensemble_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine accept_ensemble

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_accept_impl(module):accept_ensemble_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_shape_truthful_clustering_accept(module):accept_ensemble]] does both.
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
    pure subroutine accept_ensemble_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth; sizes U_history/d_history below
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_accept_impl(module):tox_stc_accept_ensemble_svd_workspace]].
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_first
            !! Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
            !! `misc/mod_STC.md`, "Output"
        integer(int32), intent(in) :: d_first
            !! Ensemble's intrinsic dimension at its bootstrap iteration
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions, n_dimensions, o), intent(in) :: U_history
            !! Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
            !! valid (see `history_len`)
        integer(int32), dimension(o), intent(in) :: d_history
            !! Trailing intrinsic dimensions, one per column of U_history; only entries
            !! 1:history_len are valid
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: history_len
            !! Number of valid columns in U_history/d_history; column history_len is the most
            !! recently accepted iteration
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `o`.
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at the most recently accepted iteration
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: normal_error_t
            !! Ensemble's normal_error at the most recently accepted iteration, see
            !! `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
            !! noise at all in its normal directions -- unlike G_t, which is a ratio already
            !! protected by its own +epsilon denominator (see `observable`) and so is
            !! required to be strictly positive.
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_tp1
            !! Candidate ensemble's tangent+normal basis
        integer(int32), intent(in) :: d_tp1
            !! Candidate ensemble's intrinsic dimension
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), intent(in) :: G_tp1
            !! Candidate ensemble's spectral gap
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: normal_error_tp1
            !! Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in) :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in) :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(out) :: tmp_m
            !! Workspace: M = U_r(:,1:d)^T U_tp1(:,1:d) for whichever reference is being compared
        real(real64), dimension(n_dimensions), intent(out) :: tmp_s
            !! Workspace: singular values of tmp_m (= cos of the principal angles)
        real(real64), dimension(lwork), intent(out) :: tmp_work
            !! Workspace: LAPACK dgesvd scratch
        logical(c_bool), intent(out) :: is_accepted
            !! .true. if all four acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=1_int32)
        call validate_in_range_int(o, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(d_first, ierr, arg_pos=4_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_int(history_len, ierr, arg_pos=7_int32, min=1_int32, max=o)
        call validate_in_range_real(G_t, ierr, arg_pos=8_int32, min=above(0.0_real64))
        call validate_in_range_real(normal_error_t, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_in_range_int(d_tp1, ierr, arg_pos=11_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_tp1, ierr, arg_pos=12_int32, min=above(0.0_real64))
        call validate_in_range_real(normal_error_tp1, ierr, arg_pos=13_int32, min=0.0_real64)
        call validate_in_range_real(chordal_dist_max_as_prcnt_of_range, ierr, arg_pos=14_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(d_max, ierr, arg_pos=15_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=16_int32, min=0.0_real64)
        call validate_in_range_real(RMSE_change_max, ierr, arg_pos=17_int32, min=0.0_real64)
        call validate_dimension_size(lwork, ierr, arg_pos=18_int32)
        call validate_all_in_range_real(U_first, n_dimensions * n_dimensions, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(U_history, n_dimensions * n_dimensions * o, ierr, arg_pos=5_int32)
        call validate_all_in_range_int(d_history, o, ierr, arg_pos=6_int32, min=0_int32, max=n_dimensions)
        call validate_all_in_range_real(U_tp1, n_dimensions * n_dimensions, ierr, arg_pos=10_int32)
        if (is_err(ierr)) return
#endif

        call accept_ensemble_impl(&
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
        call clear_err_arg_pos(ierr)
    end subroutine accept_ensemble_expert

end module tox_shape_truthful_clustering_accept
