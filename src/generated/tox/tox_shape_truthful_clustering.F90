#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering
    use tox_shape_truthful_clustering_kernel, only: ensemble_identification_kernel, ensemble_identification_merged_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math, only: above
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, clear_err_arg_pos
    use tox_errors, only: set_err_once, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: ensemble_identification
    public :: ensemble_identification_merged

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_kernel(module):ensemble_identification_kernel]].
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
    subroutine ensemble_identification(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_index,&
            k_min,&
            alpha_max,&
            d_max,&
            G_max,&
            f_max,&
            a,&
            o,&
            final_ensemble_mask,&
            stop_reason,&
            growth_radius,&
            U_history,&
            S_history,&
            d_history,&
            G_history,&
            mu_history,&
            k_history,&
            accepted_history,&
            member_added_at_step,&
            low_confidence_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth (`misc/mod_STC.md` suggests 10 as a
            !! sensible default). Always required, never optional with an auto-applied
            !! default here: a Fortran array bound cannot depend on a possibly-absent
            !! optional dummy, and this argument sizes every history output below.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to grow an ensemble around
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size for this seed's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians), see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! The minimum valid value is `above(0.0_real64)`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        logical, dimension(n_vectors), intent(out) :: final_ensemble_mask
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
        real(real64), dimension(n_dimensions, n_dimensions, o), intent(out) :: U_history
            !! Trailing tangent+normal bases, one per retained iteration, oldest to newest;
            !! zero beyond the number of iterations actually retained, see `k_history`
        real(real64), dimension(n_dimensions, o), intent(out) :: S_history
            !! Trailing singular values -- not eigenvalues, see "Output" in `misc/mod_STC.md`
            !! -- zero-padded beyond rank and beyond the number of retained iterations
        integer(int32), dimension(o), intent(out) :: d_history
            !! Trailing intrinsic dimensions, one per retained iteration
        real(real64), dimension(o), intent(out) :: G_history
            !! Trailing spectral gaps, one per retained iteration
        real(real64), dimension(n_dimensions, o), intent(out) :: mu_history
            !! Trailing ensemble centers, one per retained iteration
        integer(int32), dimension(o), intent(out) :: k_history
            !! Trailing ensemble sizes, one per retained iteration. 0 marks a column beyond
            !! the number of iterations actually retained -- a real ensemble size is always
            !! at least 1.
        logical, dimension(o), intent(out) :: accepted_history
            !! Whether the growth iteration retained in the corresponding column was
            !! accepted. Iteration 1 (the bootstrap step) is always `.true.` by convention.
            !! The single most recent column is `.false.` when, and only when, growth
            !! stopped via `STOP_REASON_REJECTED_AFTER_STABLE` or
            !! `STOP_REASON_REJECTED_IMMEDIATELY` -- see the module-level note above.
        integer(int32), dimension(n_vectors), intent(out) :: member_added_at_step
            !! `MEMBER_ADDED_AT_STEP_NON_MEMBER` for non-members, `MEMBER_ADDED_AT_STEP_SEED`
            !! for the seed itself, the growth-iteration index at which each other member
            !! joined otherwise
        logical, dimension(n_vectors), intent(out) :: low_confidence_mask
            !! Membership from this seed's iteration 1 (the unconditional bootstrap
            !! grow_ensemble+observable call), reported regardless of stop_reason -- including
            !! when stop_reason is STOP_REASON_MAX_SIZE, for which final_ensemble_mask is
            !! all-.false. All-.false. here too whenever iteration 1 itself never produced a
            !! genuine observable (an isolated seed, or a seed whose very first growth step
            !! already exceeds f_max*N) -- see "Ensemble identification", "Output" in
            !! misc/mod_STC.md
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only on a genuine LAPACK SVD non-convergence
            !! in `observable`/`accept_ensemble` -- every Stop Condition is a valid,
            !! non-error algorithmic outcome, see `misc/mod_STC.md`, "Stop Conditions".

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(seed_index, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        call validate_in_range_int(k_min, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(alpha_max, ierr, arg_pos=8_int32, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, arg_pos=9_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_in_range_real(f_max, ierr, arg_pos=11_int32, min=above(0.0_real64), max=1.0_real64)
        call validate_in_range_int(a, ierr, arg_pos=12_int32, min=1_int32)
        call validate_in_range_int(o, ierr, arg_pos=13_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call ensemble_identification_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_index = seed_index,&
            k_min = k_min,&
            alpha_max = alpha_max,&
            d_max = d_max,&
            G_max = G_max,&
            f_max = f_max,&
            a = a,&
            o = o,&
            final_ensemble_mask = final_ensemble_mask,&
            stop_reason = stop_reason,&
            growth_radius = growth_radius,&
            U_history = U_history,&
            S_history = S_history,&
            d_history = d_history,&
            G_history = G_history,&
            mu_history = mu_history,&
            k_history = k_history,&
            accepted_history = accepted_history,&
            member_added_at_step = member_added_at_step,&
            low_confidence_mask = low_confidence_mask,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine ensemble_identification

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_kernel(module):ensemble_identification_merged_kernel]].
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
    subroutine ensemble_identification_merged(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_selection_mask,&
            n_selected_seed,&
            k_min,&
            alpha_max,&
            d_max,&
            G_max,&
            f_max,&
            a,&
            o,&
            ensemble_masks,&
            ensemble_stop_reason,&
            ensemble_growth_radii,&
            ensemble_U_history,&
            ensemble_S_history,&
            ensemble_d_history,&
            ensemble_G_history,&
            ensemble_mu_history,&
            ensemble_k_history,&
            ensemble_accepted_history,&
            ensemble_member_added_at_step,&
            ensemble_low_confidence_masks,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_selected_seed
            !! Number of selected seeds (count of .TRUE. in seed_selection_mask); zero is a
            !! valid, well-defined "no ensembles" input, not an error
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth, see `ensemble_identification`. Always
            !! required, for the same reason as there: it sizes every history output below.
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        logical, dimension(n_vectors), intent(in) :: seed_selection_mask
            !! Seed selection, see `seeds`
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size for each seed's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians), see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(in), optional :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! The minimum valid value is `above(0.0_real64)`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        logical, dimension(n_vectors, n_selected_seed), intent(out) :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed, see `final_ensemble_mask`
        integer(int32), dimension(n_selected_seed), intent(out) :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `ensemble_identification`
        real(real64), dimension(n_selected_seed), intent(out) :: ensemble_growth_radii
            !! Per-ensemble growth radius, see "Local Radius Identification"
        real(real64), dimension(n_dimensions, n_dimensions, o, n_selected_seed), intent(out) :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see `U_history`
        real(real64), dimension(n_dimensions, o, n_selected_seed), intent(out) :: ensemble_S_history
            !! Per-ensemble trailing singular values, see `S_history`
        integer(int32), dimension(o, n_selected_seed), intent(out) :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions, see `d_history`
        real(real64), dimension(o, n_selected_seed), intent(out) :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps, see `G_history`
        real(real64), dimension(n_dimensions, o, n_selected_seed), intent(out) :: ensemble_mu_history
            !! Per-ensemble trailing centers, see `mu_history`
        integer(int32), dimension(o, n_selected_seed), intent(out) :: ensemble_k_history
            !! Per-ensemble trailing sizes, see `k_history`
        logical, dimension(o, n_selected_seed), intent(out) :: ensemble_accepted_history
            !! Per-ensemble trailing accepted flags, see `accepted_history`
        integer(int32), dimension(n_vectors, n_selected_seed), intent(out) :: ensemble_member_added_at_step
            !! Per-ensemble growth-iteration-joined bookkeeping, see `member_added_at_step`
        logical, dimension(n_vectors, n_selected_seed), intent(out) :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership, see `low_confidence_mask`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a genuine LAPACK SVD non-convergence
            !! occurred for any seed -- see `ensemble_identification`.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_seed, ierr, arg_pos=7_int32, min=0_int32, max=n_vectors)
        call validate_in_range_int(k_min, ierr, arg_pos=8_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(alpha_max, ierr, arg_pos=9_int32, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, arg_pos=10_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_real(f_max, ierr, arg_pos=12_int32, min=above(0.0_real64), max=1.0_real64)
        call validate_in_range_int(a, ierr, arg_pos=13_int32, min=1_int32)
        call validate_in_range_int(o, ierr, arg_pos=14_int32, min=1_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (count(seed_selection_mask, kind=int32) /= n_selected_seed) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return
#endif

        call ensemble_identification_merged_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_selection_mask = seed_selection_mask,&
            n_selected_seed = n_selected_seed,&
            k_min = k_min,&
            alpha_max = alpha_max,&
            d_max = d_max,&
            G_max = G_max,&
            f_max = f_max,&
            a = a,&
            o = o,&
            ensemble_masks = ensemble_masks,&
            ensemble_stop_reason = ensemble_stop_reason,&
            ensemble_growth_radii = ensemble_growth_radii,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history,&
            ensemble_member_added_at_step = ensemble_member_added_at_step,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine ensemble_identification_merged

end module tox_shape_truthful_clustering
