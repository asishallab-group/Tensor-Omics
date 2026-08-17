#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering(module)]]
!| # Shape Truthful Clustering (STC)
!|
!| A renormalization-group-inspired ensemble-growth clustering method; see
!| `misc/mod_STC.md` for the full algorithm definition. The family is split over five
!| kernel modules -- seeding, ensemble growing, observable, accept, and reconciliation --
!| which this module gathers, so a caller reaches the whole pipeline through it.
!|
!| Unlike `codegen_guide.md` section 5.15's own `tox_data_integration_impl` example, this
!| parent also holds `ensemble_identification`'s own kernel directly -- it is this family's
!| natural top-level entry point (it orchestrates every other SKG here), not just a bag of
!| siblings.
!|
!| `ensemble_identification_impl` grows and tracks a single ensemble from a single seed --
!| see `misc/mod_STC.md`, "Ensemble identification", "### Output". `ensemble_identification_merged_impl`
!| calls it once per seed and assembles the "#### Merged output" arrays (`ensemble_masks`,
!| `ensemble_U_history`, ...). `ensemble_reconciliation_impl`, in the sibling
!| `tox_shape_truthful_clustering_reconciliation_impl` module, then identifies and groups
!| intersecting ensembles from that merged output, on the side.
module tox_shape_truthful_clustering_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: ensemble_identification_c
    public :: ensemble_identification_merged_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering(module):ensemble_identification(subroutine)]]
    !| The iteration wrapper described in `misc/mod_STC.md`, "Ensemble identification": each
    !| growth step is `grow_ensemble` + `observable`, compared against the last *accepted*
    !| iteration via `accept_ensemble` -- skipped, by convention, for the very first growth
    !| step, see "First growth step" in the spec -- until one of the four documented Stop
    !| Conditions applies. Deliberately sequential: growth at each step depends on the
    !| previous one, so there is nothing to parallelize *within* a single seed's growth --
    !| outer-level parallelism belongs in whatever calls this once per seed, matching
    !| `grow_ensemble_impl`'s own precedent.
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
    subroutine ensemble_identification_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_index,&
            k_min,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
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
            U_first,&
            d_first,&
            ierr&
        ) bind(C, name="ensemble_identification_c")
        use tox_shape_truthful_clustering, only: ensemble_identification

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth (`misc/mod_STC.md` suggests 10 as a
            !! sensible default). Always required, never optional with an auto-applied
            !! default here: a Fortran array bound cannot depend on a possibly-absent
            !! optional dummy, and this argument sizes every history output below.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to grow an ensemble around
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: k_min
            !! Neighborhood size for this seed's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! The minimum valid value is `above(0.0_real64)`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(in), target :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        logical(c_bool), dimension(n_vectors), intent(out), target :: final_ensemble_mask
            !! The last accepted ensemble's membership. All `.false.` when `stop_reason` is
            !! `STOP_REASON_MAX_SIZE` -- see Stop Condition 1.
        integer(c_int), intent(out), target :: stop_reason
            !! Which Stop Condition ended growth: one of
            !! [[tox_shape_truthful_clustering_impl(module):STOP_REASON_MAX_SIZE(variable)]],
            !! [[tox_shape_truthful_clustering_impl(module):STOP_REASON_REJECTED_AFTER_STABLE(variable)]],
            !! [[tox_shape_truthful_clustering_impl(module):STOP_REASON_REJECTED_IMMEDIATELY(variable)]], or
            !! [[tox_shape_truthful_clustering_impl(module):STOP_REASON_FIXED_POINT(variable)]] --
            !! or [[tox_shape_truthful_clustering_impl(module):STOP_REASON_ERROR(variable)]] if
            !! `ierr` is non-zero, in which case every other output for this seed is undefined.
        real(c_double), intent(out), target :: growth_radius
            !! This seed's growth radius, see `calc_ensemble_growth_radius`
        real(c_double), dimension(n_dimensions, n_dimensions, o), intent(out), target :: U_history
            !! Trailing tangent+normal bases, one per retained iteration, oldest to newest;
            !! zero beyond the number of iterations actually retained, see `k_history`
        real(c_double), dimension(n_dimensions, o), intent(out), target :: S_history
            !! Trailing singular values -- not eigenvalues, see "Output" in `misc/mod_STC.md`
            !! -- zero-padded beyond rank and beyond the number of retained iterations
        integer(c_int), dimension(o), intent(out), target :: d_history
            !! Trailing intrinsic dimensions, one per retained iteration
        real(c_double), dimension(o), intent(out), target :: G_history
            !! Trailing spectral gaps, one per retained iteration
        real(c_double), dimension(n_dimensions, o), intent(out), target :: mu_history
            !! Trailing ensemble centers, one per retained iteration
        integer(c_int), dimension(o), intent(out), target :: k_history
            !! Trailing ensemble sizes, one per retained iteration. 0 marks a column beyond
            !! the number of iterations actually retained -- a real ensemble size is always
            !! at least 1.
        logical(c_bool), dimension(o), intent(out), target :: accepted_history
            !! Whether the growth iteration retained in the corresponding column was
            !! accepted. Iteration 1 (the bootstrap step) is always `.true.` by convention.
            !! The single most recent column is `.false.` when, and only when, growth
            !! stopped via `STOP_REASON_REJECTED_AFTER_STABLE` or
            !! `STOP_REASON_REJECTED_IMMEDIATELY` -- see the module-level note above.
        integer(c_int), dimension(n_vectors), intent(out), target :: member_added_at_step
            !! `MEMBER_ADDED_AT_STEP_NON_MEMBER` for non-members, `MEMBER_ADDED_AT_STEP_SEED`
            !! for the seed itself, the growth-iteration index at which each other member
            !! joined otherwise
        logical(c_bool), dimension(n_vectors), intent(out), target :: low_confidence_mask
            !! Membership from this seed's iteration 1 (the unconditional bootstrap
            !! grow_ensemble+observable call), reported regardless of stop_reason -- including
            !! when stop_reason is STOP_REASON_MAX_SIZE, for which final_ensemble_mask is
            !! all-.false. All-.false. here too whenever iteration 1 itself never produced a
            !! genuine observable (an isolated seed, or a seed whose very first growth step
            !! already exceeds f_max*N) -- see "Ensemble identification", "Output" in
            !! misc/mod_STC.md
        real(c_double), dimension(n_dimensions, n_dimensions), intent(out), target :: U_first
            !! Tangent+normal basis at the bootstrap iteration (iteration 1), retained for the
            !! whole growth, never evicted by the trailing o-window above -- see
            !! `accept_ensemble`'s tangent-space-drift criterion in misc/mod_STC.md. All-zero
            !! whenever iteration 1 itself never produced a genuine observable, same condition
            !! as `low_confidence_mask` above.
        integer(c_int), intent(out), target :: d_first
            !! Intrinsic dimension at the bootstrap iteration, see `U_first`. Zero under the
            !! same all-zero condition as `U_first`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success. Set only on a genuine LAPACK SVD non-convergence
            !! in `observable`/`accept_ensemble` -- every Stop Condition is a valid,
            !! non-error algorithmic outcome, see `misc/mod_STC.md`, "Stop Conditions".

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(seed_index)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(f_max)
        M_CHECK_NON_NULL(a)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(stop_reason)
        M_CHECK_NON_NULL(growth_radius)
        M_CHECK_NON_NULL(d_first)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(final_ensemble_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(U_history, n_dimensions * n_dimensions * o)
        M_CHECK_ARRAY_NON_NULL(S_history, n_dimensions * o)
        M_CHECK_ARRAY_NON_NULL(d_history, o)
        M_CHECK_ARRAY_NON_NULL(G_history, o)
        M_CHECK_ARRAY_NON_NULL(mu_history, n_dimensions * o)
        M_CHECK_ARRAY_NON_NULL(k_history, o)
        M_CHECK_ARRAY_NON_NULL(accepted_history, o)
        M_CHECK_ARRAY_NON_NULL(member_added_at_step, n_vectors)
        M_CHECK_ARRAY_NON_NULL(low_confidence_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(U_first, n_dimensions * n_dimensions)

        call ensemble_identification(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_index = seed_index,&
            k_min = k_min,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
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
            U_first = U_first,&
            d_first = d_first,&
            ierr = ierr&
        )
    end subroutine ensemble_identification_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering(module):ensemble_identification_merged(subroutine)]]
    !| See `misc/mod_STC.md`, "Ensemble identification", "#### Merged output": every per-seed
    !| output of `ensemble_identification` gains one extra trailing dimension of size
    !| `n_selected_seed`, one column per seed, in the order seeds occur in `seed_selection_mask`. Each
    !| seed's growth is fully independent of every other's, so the per-seed calls below write
    !| to disjoint array sections -- safe for `do concurrent`, matching the spec's own "In
    !| parallel grow ensembles around each seed vector" and this codebase's existing
    !| `do concurrent`-everywhere convention (`tox_shift_vectors_impl`, `tox_gene_centroids_impl`,
    !| ...). Left for a later pass if it turns out to matter in practice: the spec's own caveat
    !| that `do concurrent` is unsafe together with external-library calls under gfortran --
    !| `ensemble_identification_impl` calls LAPACK (`dgesdd`/`dgesvd`) by way of `observable`
    !| and `accept_ensemble` -- has not been stress-tested here; `!$omp parallel do` is the
    !| documented fallback if it ever is.
    !|
    !| `ensemble_member_added_at_step` is always collected (see `misc/mod_STC.md`'s "optional,
    !| user flag decides" note): unlike Ensemble Reconciliation's JSI, which the user
    !| explicitly required to be gated because it adds a real, if small, extra cost per pair,
    !| this is just bookkeeping already computed as a side effect of the per-seed growth loop
    !| itself -- there is no separate cost left to gate.
    subroutine ensemble_identification_merged_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_selection_mask,&
            n_selected_seed,&
            k_min,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
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
            ensemble_U_first,&
            ensemble_d_first,&
            ierr&
        ) bind(C, name="ensemble_identification_merged_c")
        use tox_shape_truthful_clustering, only: ensemble_identification_merged

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_seed
            !! Number of selected seeds (count of .TRUE. in seed_selection_mask); zero is a
            !! valid, well-defined "no ensembles" input, not an error
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth, see `ensemble_identification`. Always
            !! required, for the same reason as there: it sizes every history output below.
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        logical(c_bool), dimension(n_vectors), intent(in), target :: seed_selection_mask
            !! Seed selection, see `seeds`
        integer(c_int), intent(in), target :: k_min
            !! Neighborhood size for each seed's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! Maximum tolerated chordal distance between tangent bases, as a fraction of its
            !! own [0, sqrt(d)] range, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: RMSE_change_max
            !! Maximum tolerated |log(RMSE_tp1/RMSE_t)|, see `accept_ensemble`
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: f_max
            !! Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
            !! The minimum valid value is `above(0.0_real64)`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(in), target :: a
            !! Minimum accepted-iteration count for a later rejection to count as "stable", see
            !! Stop Condition 2
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(out), target :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed, see `final_ensemble_mask`
        integer(c_int), dimension(n_selected_seed), intent(out), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `ensemble_identification`
        real(c_double), dimension(n_selected_seed), intent(out), target :: ensemble_growth_radii
            !! Per-ensemble growth radius, see "Local Radius Identification"
        real(c_double), dimension(n_dimensions, n_dimensions, o, n_selected_seed), intent(out), target :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see `U_history`
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(out), target :: ensemble_S_history
            !! Per-ensemble trailing singular values, see `S_history`
        integer(c_int), dimension(o, n_selected_seed), intent(out), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions, see `d_history`
        real(c_double), dimension(o, n_selected_seed), intent(out), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps, see `G_history`
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(out), target :: ensemble_mu_history
            !! Per-ensemble trailing centers, see `mu_history`
        integer(c_int), dimension(o, n_selected_seed), intent(out), target :: ensemble_k_history
            !! Per-ensemble trailing sizes, see `k_history`
        logical(c_bool), dimension(o, n_selected_seed), intent(out), target :: ensemble_accepted_history
            !! Per-ensemble trailing accepted flags, see `accepted_history`
        integer(c_int), dimension(n_vectors, n_selected_seed), intent(out), target :: ensemble_member_added_at_step
            !! Per-ensemble growth-iteration-joined bookkeeping, see `member_added_at_step`
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(out), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership, see `low_confidence_mask`
        real(c_double), dimension(n_dimensions, n_dimensions, n_selected_seed), intent(out), target :: ensemble_U_first
            !! Per-ensemble bootstrap-iteration tangent+normal basis, see `U_first`
        integer(c_int), dimension(n_selected_seed), intent(out), target :: ensemble_d_first
            !! Per-ensemble bootstrap-iteration intrinsic dimension, see `d_first`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success. Set only if a genuine LAPACK SVD non-convergence
            !! occurred for any seed -- see `ensemble_identification`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_seed)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(f_max)
        M_CHECK_NON_NULL(a)
        M_CHECK_NON_NULL(o)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(seed_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_growth_radii, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_history, n_dimensions * n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_accepted_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_member_added_at_step, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_first, n_dimensions * n_dimensions * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_first, n_selected_seed)

        call ensemble_identification_merged(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_selection_mask = seed_selection_mask,&
            n_selected_seed = n_selected_seed,&
            k_min = k_min,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
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
            ensemble_U_first = ensemble_U_first,&
            ensemble_d_first = ensemble_d_first,&
            ierr = ierr&
        )
    end subroutine ensemble_identification_merged_c

end module tox_shape_truthful_clustering_c
#endif
