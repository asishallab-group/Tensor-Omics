#include <src/macros.h>

!> # Shape Truthful Clustering (STC)
!|
!| A renormalization-group-inspired ensemble-growth clustering method; see
!| `misc/mod_STC.md` for the full algorithm definition. The family is split over five
!| kernel modules -- seeding, ensemble growing, observable, accept, and reconciliation --
!| which this module gathers, so a caller reaches the whole pipeline through it.
!|
!| Unlike `codegen_guide.md` section 5.15's own `tox_data_integration_kernel` example, this
!| parent also holds `ensemble_identification`'s own kernel directly -- it is this family's
!| natural top-level entry point (it orchestrates every other SKG here), not just a bag of
!| siblings.
!|
!| `ensemble_identification_kernel` grows and tracks a single ensemble from a single seed --
!| see `misc/mod_STC.md`, "Ensemble identification", "### Output". `ensemble_identification_merged_kernel`
!| calls it once per seed and assembles the "#### Merged output" arrays (`ensemble_masks`,
!| `ensemble_U_history`, ...). `ensemble_reconciliation_kernel`, in the sibling
!| `tox_shape_truthful_clustering_reconciliation_kernel` module, then identifies and groups
!| intersecting ensembles from that merged output, on the side.
module tox_shape_truthful_clustering_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err_once, ERR_INTERNAL
    use f42_utils, only: above
    use tox_shape_truthful_clustering_seeding_kernel
    use tox_shape_truthful_clustering_ensemble_growing_kernel
    use tox_shape_truthful_clustering_observable_kernel
    use tox_shape_truthful_clustering_accept_kernel
    use tox_shape_truthful_clustering_reconciliation_kernel
    M_IMPLICIT_NONE

#define CM_STC_GROWTH_RADIUS_K_MIN_DEFAULT 30_int32
#define CM_STC_F_MAX_DEFAULT 0.95_real64
#define CM_STC_A_DEFAULT 2_int32
#define CM_STC_STOP_REASON_MAX_SIZE 1_int32
#define CM_STC_STOP_REASON_REJECTED_AFTER_STABLE 2_int32
#define CM_STC_STOP_REASON_REJECTED_IMMEDIATELY 3_int32
#define CM_STC_STOP_REASON_FIXED_POINT 4_int32
#define CM_STC_STOP_REASON_ERROR 0_int32
#define CM_STC_MEMBER_ADDED_AT_STEP_NON_MEMBER (-1_int32)
#define CM_STC_MEMBER_ADDED_AT_STEP_SEED 0_int32

    private

    !> Stop Condition 1, see `misc/mod_STC.md`, "Stop Conditions": maximum ensemble size
    !| reached -- no ensemble is returned for the seed.
    integer(int32), parameter, public :: STOP_REASON_MAX_SIZE = CM_STC_STOP_REASON_MAX_SIZE
    !> Stop Condition 2: rejected after being stably accepted at least `a` times.
    integer(int32), parameter, public :: STOP_REASON_REJECTED_AFTER_STABLE = CM_STC_STOP_REASON_REJECTED_AFTER_STABLE
    !> Stop Condition 3: rejected on the very first `accept_ensemble` check.
    integer(int32), parameter, public :: STOP_REASON_REJECTED_IMMEDIATELY = CM_STC_STOP_REASON_REJECTED_IMMEDIATELY
    !> Stop Condition 4: `grow_ensemble` returned an unchanged membership, a natural fixed point.
    integer(int32), parameter, public :: STOP_REASON_FIXED_POINT = CM_STC_STOP_REASON_FIXED_POINT
    !> Not a Stop Condition: growth did not complete, because a genuine LAPACK SVD
    !| non-convergence set `ierr`. This is the value `stop_reason` holds for a seed until one
    !| of the four Stop Conditions above is actually reached -- so a caller that reads
    !| `stop_reason` without checking `ierr` first cannot mistake an aborted seed for
    !| Stop Condition 4's own, legitimate fixed point (both would otherwise default to the
    !| same value).
    integer(int32), parameter, public :: STOP_REASON_ERROR = CM_STC_STOP_REASON_ERROR
    !> `member_added_at_step` sentinel for a vector that never joined the final ensemble.
    integer(int32), parameter, public :: MEMBER_ADDED_AT_STEP_NON_MEMBER = CM_STC_MEMBER_ADDED_AT_STEP_NON_MEMBER
    !> `member_added_at_step` value for the seed vector itself (joined "at step 0", before growth).
    integer(int32), parameter, public :: MEMBER_ADDED_AT_STEP_SEED = CM_STC_MEMBER_ADDED_AT_STEP_SEED

    public :: ensemble_identification_kernel
    public :: ensemble_identification_merged_kernel

contains

    !> summary: Grow and track a single ensemble from one seed until a Stop Condition is reached
    !| AUTHOR_ASIS_HALLAB
    !| The iteration wrapper described in `misc/mod_STC.md`, "Ensemble identification": each
    !| growth step is `grow_ensemble` + `observable`, compared against the last *accepted*
    !| iteration via `accept_ensemble` -- skipped, by convention, for the very first growth
    !| step, see "First growth step" in the spec -- until one of the four documented Stop
    !| Conditions applies. Deliberately sequential: growth at each step depends on the
    !| previous one, so there is nothing to parallelize *within* a single seed's growth --
    !| outer-level parallelism belongs in whatever calls this once per seed, matching
    !| `grow_ensemble_kernel`'s own precedent.
    !|
    !| Two deliberate readings of the spec, flagged here since the prose leaves them
    !| implicit: (1) an isolated seed -- no neighbor at all within its own growth radius --
    !| is reported as Stop Condition 4 (a natural fixed point) with a trivial one-member
    !| `final_ensemble_mask`, not as "no ensemble" (that is Stop Condition 1's own, distinct
    !| meaning). (2) On a rejection (Stop Condition 2 or 3), the *rejected* candidate's
    !| observable is still pushed into the trailing history (marked `.false.` in
    !| `accepted_history`) so a caller can see what got rejected and why, even though
    !| `final_ensemble_mask` reflects the last *accepted* state, not this one -- otherwise
    !| `accepted_history` could only ever read `.true.`, since only accepted iterations would
    !| ever reach the array at all.
    pure subroutine ensemble_identification_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                    seed_index, k_min, alpha_max, d_max, G_max, f_max, a, o, &
                                                    final_ensemble_mask, stop_reason, growth_radius, &
                                                    U_history, S_history, d_history, G_history, mu_history, &
                                                    k_history, accepted_history, member_added_at_step, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree index over `vectors`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Dimension order used to build `kd_indices`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to grow an ensemble around
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size for this seed's growth radius, see `calc_ensemble_growth_radius`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
            !! DM_DEFAULT(CM_STC_GROWTH_RADIUS_K_MIN_DEFAULT)
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians), see `accept_ensemble`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(2.0_real64 * atan(1.0_real64))
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! DM_MIN(above(0.0_real64))
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_STC_F_MAX_DEFAULT)
        integer(int32), intent(in), optional :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_STC_A_DEFAULT)
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth (`misc/mod_STC.md` suggests 10 as a
            !! sensible default). Always required, never optional with an auto-applied
            !! default here: a Fortran array bound cannot depend on a possibly-absent
            !! optional dummy, and this argument sizes every history output below.
            !! DM_MIN(1_int32)
        logical, intent(out) :: final_ensemble_mask(n_vectors)
            !! The last accepted ensemble's membership. All `.false.` when `stop_reason` is
            !! `STOP_REASON_MAX_SIZE` -- see Stop Condition 1.
        integer(int32), intent(out) :: stop_reason
            !! Which Stop Condition ended growth: one of
            !! [[tox_shape_truthful_clustering_kernel(module):STOP_REASON_MAX_SIZE(variable)]],
            !! [[tox_shape_truthful_clustering_kernel(module):STOP_REASON_REJECTED_AFTER_STABLE(variable)]],
            !! [[tox_shape_truthful_clustering_kernel(module):STOP_REASON_REJECTED_IMMEDIATELY(variable)]], or
            !! [[tox_shape_truthful_clustering_kernel(module):STOP_REASON_FIXED_POINT(variable)]] --
            !! or [[tox_shape_truthful_clustering_kernel(module):STOP_REASON_ERROR(variable)]] if
            !! `ierr` is non-zero, in which case every other output for this seed is undefined.
        real(real64), intent(out) :: growth_radius
            !! This seed's growth radius, see `calc_ensemble_growth_radius`
        real(real64), intent(out) :: U_history(n_dimensions, n_dimensions, o)
            !! Trailing tangent+normal bases, one per retained iteration, oldest to newest;
            !! zero beyond the number of iterations actually retained, see `k_history`
        real(real64), intent(out) :: S_history(n_dimensions, o)
            !! Trailing singular values -- not eigenvalues, see "Output" in `misc/mod_STC.md`
            !! -- zero-padded beyond rank and beyond the number of retained iterations
        integer(int32), intent(out) :: d_history(o)
            !! Trailing intrinsic dimensions, one per retained iteration
        real(real64), intent(out) :: G_history(o)
            !! Trailing spectral gaps, one per retained iteration
        real(real64), intent(out) :: mu_history(n_dimensions, o)
            !! Trailing ensemble centers, one per retained iteration
        integer(int32), intent(out) :: k_history(o)
            !! Trailing ensemble sizes, one per retained iteration. 0 marks a column beyond
            !! the number of iterations actually retained -- a real ensemble size is always
            !! at least 1.
        logical, intent(out) :: accepted_history(o)
            !! Whether the growth iteration retained in the corresponding column was
            !! accepted. Iteration 1 (the bootstrap step) is always `.true.` by convention.
            !! The single most recent column is `.false.` when, and only when, growth
            !! stopped via `STOP_REASON_REJECTED_AFTER_STABLE` or
            !! `STOP_REASON_REJECTED_IMMEDIATELY` -- see the module-level note above.
        integer(int32), intent(out) :: member_added_at_step(n_vectors)
            !! `MEMBER_ADDED_AT_STEP_NON_MEMBER` for non-members, `MEMBER_ADDED_AT_STEP_SEED`
            !! for the seed itself, the growth-iteration index at which each other member
            !! joined otherwise
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only on a genuine LAPACK SVD non-convergence
            !! in `observable`/`accept_ensemble` -- every Stop Condition is a valid,
            !! non-error algorithmic outcome, see `misc/mod_STC.md`, "Stop Conditions".

        integer(int32) :: growth_neighbors(n_vectors)
        real(real64)    :: growth_distances(n_vectors)
        integer(int32) :: range_stack(3, n_vectors)
        integer(int32) :: growth_sort_perm(n_vectors)
        logical         :: member_mask_buf(n_vectors)

        logical         :: mask_current(n_vectors)
        logical         :: mask_candidate(n_vectors)
        integer(int32) :: n_candidate

        real(real64)   :: U_candidate(n_dimensions, n_dimensions)
        real(real64)   :: eigenvalues_candidate(n_dimensions)
        real(real64)   :: mu_candidate(n_dimensions)
        real(real64)   :: S_candidate(n_dimensions)
        real(real64)   :: tangent_scales_scratch(n_dimensions)
        real(real64)   :: normal_error_scratch
        integer(int32) :: d_candidate
        real(real64)   :: G_candidate
        integer(int32) :: obs_ierr

        real(real64)   :: prev_U(n_dimensions, n_dimensions)
        integer(int32) :: prev_d
        real(real64)   :: prev_G

        logical         :: is_accepted
        integer(int32) :: accept_ierr

        real(real64)   :: actual_f_max
        integer(int32) :: actual_a
        integer(int32) :: accepted_count, history_len, t

        call set_ok(ierr)

        M_DEFAULT_VAL(f_max, actual_f_max, CM_STC_F_MAX_DEFAULT)
        M_DEFAULT_VAL(a, actual_a, CM_STC_A_DEFAULT)

        final_ensemble_mask   = .false.
        stop_reason           = STOP_REASON_ERROR
        U_history             = 0.0_real64
        S_history              = 0.0_real64
        d_history              = 0
        G_history              = 0.0_real64
        mu_history             = 0.0_real64
        k_history               = 0
        accepted_history        = .false.
        member_added_at_step    = MEMBER_ADDED_AT_STEP_NON_MEMBER
        member_added_at_step(seed_index) = MEMBER_ADDED_AT_STEP_SEED
        history_len   = 0
        accepted_count = 0

        call calc_ensemble_growth_radius_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                seed_index, k_min, &
                                                growth_neighbors, growth_distances, range_stack, growth_sort_perm, &
                                                growth_radius)

        mask_current = .false.
        mask_current(seed_index) = .true.

        ! First growth step (the bootstrap): unconditional, no accept_ensemble check yet --
        ! see "First growth step" in `misc/mod_STC.md`.
        call grow_ensemble_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                  mask_current, growth_radius, &
                                  range_stack, member_mask_buf, mask_candidate)

        n_candidate = count(mask_candidate)

        if (n_candidate < 2) then
            ! Isolated seed: no neighbor lies within its own growth radius, so growth can
            ! never begin. A well-defined natural fixed point (Stop Condition 4), not Stop
            ! Condition 1's "no ensemble" -- this genuinely is one, just one that cannot
            ! grow. `observable` requires >= 2 members, so no history entry is possible.
            final_ensemble_mask = mask_candidate
            stop_reason          = STOP_REASON_FIXED_POINT
            return
        end if

        if (real(n_candidate, real64) >= actual_f_max*real(n_vectors, real64)) then
            stop_reason          = STOP_REASON_MAX_SIZE
            member_added_at_step = MEMBER_ADDED_AT_STEP_NON_MEMBER
            return
        end if

        block
            real(real64)   :: y_c(n_dimensions, n_candidate)
            real(real64)   :: s_c(min(n_dimensions, n_candidate))
            real(real64)   :: u_econ_c(n_dimensions, min(n_dimensions, n_candidate))
            real(real64)   :: vt_econ_c(min(n_dimensions, n_candidate), n_candidate)
            real(real64)   :: work_c(4*min(n_dimensions, n_candidate)**2 + 7*min(n_dimensions, n_candidate))
            integer(int32) :: iwork_c(8*min(n_dimensions, n_candidate))

            call observable_kernel(vectors, n_dimensions, n_vectors, mask_candidate, n_candidate, &
                                   4*min(n_dimensions, n_candidate)**2 + 7*min(n_dimensions, n_candidate), &
                                   8*min(n_dimensions, n_candidate), &
                                   y_c, s_c, u_econ_c, vt_econ_c, work_c, iwork_c, &
                                   U_candidate, eigenvalues_candidate, mu_candidate, d_candidate, G_candidate, &
                                   normal_error_scratch, tangent_scales_scratch, obs_ierr)
        end block

        if (obs_ierr /= 0) then
            call set_err_once(ierr, ERR_INTERNAL)
            return
        end if

        S_candidate = sqrt(eigenvalues_candidate*real(n_candidate - 1, real64))

        where (mask_candidate .and. .not. mask_current) member_added_at_step = 1
        mask_current        = mask_candidate
        final_ensemble_mask = mask_current
        accepted_count       = 1

        call stc_push_ensemble_history(n_dimensions, o, U_candidate, S_candidate, d_candidate, G_candidate, &
                                       mu_candidate, n_candidate, .true., history_len, &
                                       U_history, S_history, d_history, G_history, mu_history, k_history, &
                                       accepted_history)

        prev_U = U_candidate
        prev_d = d_candidate
        prev_G = G_candidate

        growth_loop: do t = 2, n_vectors

            call grow_ensemble_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                      mask_current, growth_radius, &
                                      range_stack, member_mask_buf, mask_candidate)

            if (all(mask_candidate .eqv. mask_current)) then
                stop_reason = STOP_REASON_FIXED_POINT
                exit growth_loop
            end if

            n_candidate = count(mask_candidate)

            if (real(n_candidate, real64) >= actual_f_max*real(n_vectors, real64)) then
                ! Stop Condition 1 poisons the whole seed, not just this step: no ensemble is
                ! returned at all, even though earlier iterations had already been accepted.
                stop_reason           = STOP_REASON_MAX_SIZE
                final_ensemble_mask   = .false.
                U_history             = 0.0_real64
                S_history             = 0.0_real64
                d_history             = 0
                G_history             = 0.0_real64
                mu_history            = 0.0_real64
                k_history             = 0
                accepted_history      = .false.
                member_added_at_step  = MEMBER_ADDED_AT_STEP_NON_MEMBER
                exit growth_loop
            end if

            block
                real(real64)   :: y_c(n_dimensions, n_candidate)
                real(real64)   :: s_c(min(n_dimensions, n_candidate))
                real(real64)   :: u_econ_c(n_dimensions, min(n_dimensions, n_candidate))
                real(real64)   :: vt_econ_c(min(n_dimensions, n_candidate), n_candidate)
                real(real64)   :: work_c(4*min(n_dimensions, n_candidate)**2 + 7*min(n_dimensions, n_candidate))
                integer(int32) :: iwork_c(8*min(n_dimensions, n_candidate))

                call observable_kernel(vectors, n_dimensions, n_vectors, mask_candidate, n_candidate, &
                                       4*min(n_dimensions, n_candidate)**2 + 7*min(n_dimensions, n_candidate), &
                                       8*min(n_dimensions, n_candidate), &
                                       y_c, s_c, u_econ_c, vt_econ_c, work_c, iwork_c, &
                                       U_candidate, eigenvalues_candidate, mu_candidate, d_candidate, G_candidate, &
                                       normal_error_scratch, tangent_scales_scratch, obs_ierr)
            end block

            if (obs_ierr /= 0) then
                call set_err_once(ierr, ERR_INTERNAL)
                return
            end if

            S_candidate = sqrt(eigenvalues_candidate*real(n_candidate - 1, real64))

            block
                integer(int32) :: d_common
                real(real64)   :: m_c(min(prev_d, d_candidate), min(prev_d, d_candidate))
                real(real64)   :: s_ang_c(min(prev_d, d_candidate))
                real(real64)   :: accept_work_c(max(1_int32, 5_int32*min(prev_d, d_candidate)))

                d_common = min(prev_d, d_candidate)
                call accept_ensemble_kernel(n_dimensions, prev_U, prev_d, prev_G, U_candidate, d_candidate, G_candidate, &
                                            alpha_max, d_max, G_max, max(1_int32, 5_int32*d_common), &
                                            m_c, s_ang_c, accept_work_c, is_accepted, accept_ierr)
            end block

            if (accept_ierr /= 0) then
                call set_err_once(ierr, ERR_INTERNAL)
                return
            end if

            if (is_accepted) then
                where (mask_candidate .and. .not. mask_current) member_added_at_step = t
                mask_current        = mask_candidate
                final_ensemble_mask = mask_current
                accepted_count       = accepted_count + 1

                call stc_push_ensemble_history(n_dimensions, o, U_candidate, S_candidate, d_candidate, G_candidate, &
                                               mu_candidate, n_candidate, .true., history_len, &
                                               U_history, S_history, d_history, G_history, mu_history, k_history, &
                                               accepted_history)

                prev_U = U_candidate
                prev_d = d_candidate
                prev_G = G_candidate
            else
                ! Rejected: push the rejected candidate for diagnosis (see the module-level
                ! note above), then stop -- final_ensemble_mask is already the last accepted
                ! state, no further update needed.
                call stc_push_ensemble_history(n_dimensions, o, U_candidate, S_candidate, d_candidate, G_candidate, &
                                               mu_candidate, n_candidate, .false., history_len, &
                                               U_history, S_history, d_history, G_history, mu_history, k_history, &
                                               accepted_history)

                if (accepted_count >= actual_a) then
                    stop_reason = STOP_REASON_REJECTED_AFTER_STABLE
                else
                    stop_reason = STOP_REASON_REJECTED_IMMEDIATELY
                end if
                exit growth_loop
            end if

        end do growth_loop

    end subroutine ensemble_identification_kernel

    !> summary: Run ensemble_identification once per seed and assemble the merged, per-ensemble output arrays
    !| AUTHOR_ASIS_HALLAB
    !| See `misc/mod_STC.md`, "Ensemble identification", "#### Merged output": every per-seed
    !| output of `ensemble_identification` gains one extra trailing dimension of size
    !| `n_selected_seed`, one column per seed, in the order seeds occur in `seed_selection_mask`. Each
    !| seed's growth is fully independent of every other's, so the per-seed calls below write
    !| to disjoint array sections -- safe for `do concurrent`, matching the spec's own "In
    !| parallel grow ensembles around each seed vector" and this codebase's existing
    !| `do concurrent`-everywhere convention (`tox_shift_vectors_kernel`, `tox_gene_centroids_kernel`,
    !| ...). Left for a later pass if it turns out to matter in practice: the spec's own caveat
    !| that `do concurrent` is unsafe together with external-library calls under gfortran --
    !| `ensemble_identification_kernel` calls LAPACK (`dgesdd`/`dgesvd`) by way of `observable`
    !| and `accept_ensemble` -- has not been stress-tested here; `!$omp parallel do` is the
    !| documented fallback if it ever is.
    !|
    !| `ensemble_member_added_at_step` is always collected (see `misc/mod_STC.md`'s "optional,
    !| user flag decides" note): unlike Ensemble Reconciliation's JSI, which the user
    !| explicitly required to be gated because it adds a real, if small, extra cost per pair,
    !| this is just bookkeeping already computed as a side effect of the per-seed growth loop
    !| itself -- there is no separate cost left to gate.
    pure subroutine ensemble_identification_merged_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                           seed_selection_mask, n_selected_seed, &
                                                           k_min, alpha_max, d_max, G_max, f_max, a, o, &
                                                           ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
                                                           ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                                                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                                                           ensemble_accepted_history, ensemble_member_added_at_step, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree index over `vectors`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Dimension order used to build `kd_indices`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        logical, intent(in) :: seed_selection_mask(n_vectors)
            !! Seed selection, see `seeds`
        integer(int32), intent(in) :: n_selected_seed
            !! Number of selected seeds (count of .TRUE. in seed_selection_mask); zero is a
            !! valid, well-defined "no ensembles" input, not an error
            !! DM_MIN(0_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size for each seed's growth radius, see `calc_ensemble_growth_radius`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
            !! DM_DEFAULT(CM_STC_GROWTH_RADIUS_K_MIN_DEFAULT)
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians), see `accept_ensemble`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(2.0_real64 * atan(1.0_real64))
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! DM_MIN(above(0.0_real64))
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_STC_F_MAX_DEFAULT)
        integer(int32), intent(in), optional :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_STC_A_DEFAULT)
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth, see `ensemble_identification`. Always
            !! required, for the same reason as there: it sizes every history output below.
            !! DM_MIN(1_int32)
        logical, intent(out) :: ensemble_masks(n_vectors, n_selected_seed)
            !! Per-ensemble accepted membership, one column per seed, see `final_ensemble_mask`
        integer(int32), intent(out) :: ensemble_stop_reason(n_selected_seed)
            !! Per-ensemble Stop Condition, see `ensemble_identification`
        real(real64), intent(out) :: ensemble_growth_radii(n_selected_seed)
            !! Per-ensemble growth radius, see "Local Radius Identification"
        real(real64), intent(out) :: ensemble_U_history(n_dimensions, n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing tangent+normal bases, see `U_history`
        real(real64), intent(out) :: ensemble_S_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing singular values, see `S_history`
        integer(int32), intent(out) :: ensemble_d_history(o, n_selected_seed)
            !! Per-ensemble trailing intrinsic dimensions, see `d_history`
        real(real64), intent(out) :: ensemble_G_history(o, n_selected_seed)
            !! Per-ensemble trailing spectral gaps, see `G_history`
        real(real64), intent(out) :: ensemble_mu_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing centers, see `mu_history`
        integer(int32), intent(out) :: ensemble_k_history(o, n_selected_seed)
            !! Per-ensemble trailing sizes, see `k_history`
        logical, intent(out) :: ensemble_accepted_history(o, n_selected_seed)
            !! Per-ensemble trailing accepted flags, see `accepted_history`
        integer(int32), intent(out) :: ensemble_member_added_at_step(n_vectors, n_selected_seed)
            !! Per-ensemble growth-iteration-joined bookkeeping, see `member_added_at_step`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a genuine LAPACK SVD non-convergence
            !! occurred for any seed -- see `ensemble_identification`.

        integer(int32) :: seed_indices(n_selected_seed)
        integer(int32) :: ierr_per_seed(n_selected_seed)
        integer(int32) :: i_vec, i_e, next_slot

        call set_ok(ierr)

        next_slot = 0
        do i_vec = 1, n_vectors
            if (seed_selection_mask(i_vec)) then
                next_slot = next_slot + 1
                seed_indices(next_slot) = i_vec
            end if
        end do

        ensemble_masks               = .false.
        ensemble_stop_reason         = STOP_REASON_ERROR
        ensemble_growth_radii        = 0.0_real64
        ensemble_U_history           = 0.0_real64
        ensemble_S_history           = 0.0_real64
        ensemble_d_history           = 0
        ensemble_G_history           = 0.0_real64
        ensemble_mu_history          = 0.0_real64
        ensemble_k_history           = 0
        ensemble_accepted_history    = .false.
        ensemble_member_added_at_step = MEMBER_ADDED_AT_STEP_NON_MEMBER
        ierr_per_seed                = 0

        do concurrent (i_e=1:n_selected_seed) shared(vectors, kd_indices, dimension_order, seed_indices, &
                                                      ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
                                                      ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                                                      ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                                                      ensemble_accepted_history, ensemble_member_added_at_step, &
                                                      ierr_per_seed)
            call ensemble_identification_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                seed_indices(i_e), k_min, alpha_max, d_max, G_max, f_max, a, o, &
                                                ensemble_masks(:, i_e), ensemble_stop_reason(i_e), &
                                                ensemble_growth_radii(i_e), ensemble_U_history(:, :, :, i_e), &
                                                ensemble_S_history(:, :, i_e), ensemble_d_history(:, i_e), &
                                                ensemble_G_history(:, i_e), ensemble_mu_history(:, :, i_e), &
                                                ensemble_k_history(:, i_e), ensemble_accepted_history(:, i_e), &
                                                ensemble_member_added_at_step(:, i_e), ierr_per_seed(i_e))
        end do

        do i_e = 1, n_selected_seed
            if (ierr_per_seed(i_e) /= 0) then
                call set_err_once(ierr, ERR_INTERNAL)
                exit
            end if
        end do

    end subroutine ensemble_identification_merged_kernel

    !> Shift-and-append one observable into a trailing o-column history window, oldest to
    !| newest. Not itself a kernel (private, no wrapper): pure array bookkeeping shared by
    !| every push site in `ensemble_identification_kernel`, not part of the public contract.
    pure subroutine stc_push_ensemble_history(n_dimensions, o, U_candidate, S_candidate, d_candidate, G_candidate, &
                                              mu_candidate, n_candidate, is_accepted, history_len, &
                                              U_history, S_history, d_history, G_history, mu_history, k_history, &
                                              accepted_history)
        integer(int32), intent(in) :: n_dimensions, o, d_candidate, n_candidate
        real(real64), intent(in) :: U_candidate(n_dimensions, n_dimensions)
        real(real64), intent(in) :: S_candidate(n_dimensions)
        real(real64), intent(in) :: G_candidate
        real(real64), intent(in) :: mu_candidate(n_dimensions)
        logical, intent(in) :: is_accepted
        integer(int32), intent(inout) :: history_len
        real(real64), intent(inout) :: U_history(n_dimensions, n_dimensions, o)
        real(real64), intent(inout) :: S_history(n_dimensions, o)
        integer(int32), intent(inout) :: d_history(o)
        real(real64), intent(inout) :: G_history(o)
        real(real64), intent(inout) :: mu_history(n_dimensions, o)
        integer(int32), intent(inout) :: k_history(o)
        logical, intent(inout) :: accepted_history(o)

        if (history_len < o) then
            history_len = history_len + 1
        else
            U_history(:, :, 1:o - 1) = U_history(:, :, 2:o)
            S_history(:, 1:o - 1)    = S_history(:, 2:o)
            d_history(1:o - 1)        = d_history(2:o)
            G_history(1:o - 1)        = G_history(2:o)
            mu_history(:, 1:o - 1)    = mu_history(:, 2:o)
            k_history(1:o - 1)        = k_history(2:o)
            accepted_history(1:o - 1) = accepted_history(2:o)
        end if

        U_history(:, :, history_len) = U_candidate
        S_history(:, history_len)    = S_candidate
        d_history(history_len)        = d_candidate
        G_history(history_len)        = G_candidate
        mu_history(:, history_len)    = mu_candidate
        k_history(history_len)        = n_candidate
        accepted_history(history_len) = is_accepted

    end subroutine stc_push_ensemble_history

end module tox_shape_truthful_clustering_kernel
