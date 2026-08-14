#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_reconciliation(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_reconciliation_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: ensemble_reconciliation_c
    public :: merge_to_super_ensembles_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_reconciliation(module):ensemble_reconciliation(subroutine)]]
    !| See this module's own header comment for why this is a two-call orchestrator, and
    !| `tox_shape_truthful_clustering_filter_kernel`'s own `filter_ensembles_kernel` for the
    !| full eligibility-filtering algorithm (Stop Condition/final dimension/final variance
    !| explained, each independently optional, combined by logical AND). `ierr` is set only if
    !| `merge_to_super_ensembles_kernel` discovers a component larger than `max_group_size` --
    !| see that kernel's own doc comment.
    subroutine ensemble_reconciliation_c(&
            ensemble_masks,&
            ensemble_stop_reason,&
            n_dimensions,&
            n_vectors,&
            n_ensembles,&
            ensemble_U_history,&
            ensemble_d_history,&
            ensemble_S_history,&
            ensemble_mu_history,&
            ensemble_G_history,&
            ensemble_k_history,&
            ensemble_accepted_history,&
            o,&
            mode,&
            min_overlap_coefficient,&
            report_overlap_coefficient,&
            allowed_stop_reasons,&
            d_min,&
            d_max,&
            var_explained_min,&
            max_group_size,&
            super_ensembles,&
            n_super_ensembles,&
            super_ensembles_overlap_coefficient,&
            eligible,&
            eligible_by_stop_condition,&
            eligible_by_dimension,&
            eligible_by_var_explained,&
            ierr&
        ) bind(C, name="ensemble_reconciliation_c")
        use tox_shape_truthful_clustering_reconciliation, only: ensemble_reconciliation
        use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E, see Ensemble Identification's merged `ensemble_masks`.
            !! At least 2: with fewer, no pair can ever intersect, so there is nothing this
            !! module could report or group.
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
            !! can hold; sizes its row dimension. `misc/mod_STC.md` suggests
            !! $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
            !! optional with an auto-applied default here, for the same reason as
            !! `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
            !! possibly-absent optional dummy, and a runtime-dependent value like
            !! $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
            !! default would need to be either.
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_ensembles`.
        logical(c_bool), dimension(n_vectors, n_ensembles), intent(in), target :: ensemble_masks
            !! Per-ensemble membership, see Ensemble Identification's merged output
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `filter_ensembles_kernel`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `4_int32`.
        real(c_double), dimension(n_dimensions, n_dimensions, o, n_ensembles), intent(in), target :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see Ensemble Identification's merged
            !! output
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_S_history
            !! Per-ensemble trailing singular values
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_mu_history
            !! Per-ensemble trailing centers
        real(c_double), dimension(o, n_ensembles), intent(in), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_k_history
            !! Per-ensemble trailing sizes
        logical(c_bool), dimension(o, n_ensembles), intent(in), target :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself accepted
        character(len=1, kind=c_char), dimension(25), intent(in), target :: mode
            !! How intersections are processed
            !!
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
            !! The default value is `1_int32`.
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), intent(in), target :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see `merge_to_super_ensembles_kernel`'s own note on this being guarded, not
            !! unconditional
            !! The default value is `.false.`.
        logical(c_bool), dimension(4), intent(in), optional :: allowed_stop_reasons
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_stop_condition_kernel`
        integer(c_int), intent(in), optional :: d_min
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_dimension_kernel`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), optional :: d_max
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_dimension_kernel`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), optional :: var_explained_min
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_var_explained_kernel`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(max_group_size, n_ensembles*(n_ensembles-1)), intent(out), target :: super_ensembles
            !! One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            !! belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            !! the group's actual size, and 0 in every row of an unused trailing column beyond
            !! `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]'s
            !! own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            !! intersects) -- deliberately not divided by 2: the generator translates this
            !! specification expression close to verbatim into the Python/R bindings, where
            !! `/` on two integers is true division, not Fortran's own truncating integer
            !! division, so a literal `/2` here breaks the generated Python binding (a `float`
            !! where `np.empty`'s shape wants an `int`); see `misc/code_gen_footgun.md`. A
            !! safe, if looser, upper bound for modes 2 and 3 too, whose groups can never
            !! outnumber mode 1's own worst case.
        integer(c_int), intent(out), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            !! actually filled
        real(c_double), dimension(max_group_size-1, n_ensembles*(n_ensembles-1)), intent(out), target :: super_ensembles_overlap_coefficient
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible
            !! Combined per-ensemble eligibility actually used for merging above -- see
            !! `filter_ensembles_kernel`. Ineligible ensembles are otherwise untouched: they
            !! are never removed from `ensemble_masks` or anything else this whole family
            !! reports, only excluded from contributing a pair here.
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_stop_condition
            !! See `filter_ensembles_by_stop_condition_kernel`
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_dimension
            !! See `filter_ensembles_by_dimension_kernel`
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_var_explained
            !! See `filter_ensembles_by_var_explained_kernel`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.
        logical, dimension(n_vectors, n_ensembles) :: ensemble_masks_f
        logical, dimension(o, n_ensembles) :: ensemble_accepted_history_f
        integer(int32) :: mode_mode_f
        logical :: report_overlap_coefficient_f
        logical, dimension(:), allocatable :: allowed_stop_reasons_f
        logical, dimension(n_ensembles) :: eligible_f
        logical, dimension(n_ensembles) :: eligible_by_stop_condition_f
        logical, dimension(n_ensembles) :: eligible_by_dimension_f
        logical, dimension(n_ensembles) :: eligible_by_var_explained_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_NON_NULL(report_overlap_coefficient)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_history, n_dimensions * n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_history, n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_history, n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_accepted_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(mode, 25)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_ensembles*(n_ensembles-1)))
        M_CHECK_ARRAY_NON_NULL(super_ensembles_overlap_coefficient, (max_group_size-1) * (n_ensembles*(n_ensembles-1)))
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_stop_condition, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_dimension, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_var_explained, n_ensembles)

        ensemble_masks_f = ensemble_masks
        ensemble_accepted_history_f = ensemble_accepted_history
        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("report")
                    mode_mode_f = MODE_REPORT
                case ("merge_overlap_coefficient")
                    mode_mode_f = MODE_MERGE_OVERLAP_COEFFICIENT
                case ("merge_any")
                    mode_mode_f = MODE_MERGE_ANY
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        report_overlap_coefficient_f = report_overlap_coefficient
        if (present(allowed_stop_reasons)) allowed_stop_reasons_f = allowed_stop_reasons

        call ensemble_reconciliation(&
            ensemble_masks = ensemble_masks_f,&
            ensemble_stop_reason = ensemble_stop_reason,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            n_ensembles = n_ensembles,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history_f,&
            o = o,&
            mode = mode_mode_f,&
            min_overlap_coefficient = min_overlap_coefficient,&
            report_overlap_coefficient = report_overlap_coefficient_f,&
            allowed_stop_reasons = allowed_stop_reasons_f,&
            d_min = d_min,&
            d_max = d_max,&
            var_explained_min = var_explained_min,&
            max_group_size = max_group_size,&
            super_ensembles = super_ensembles,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles_overlap_coefficient = super_ensembles_overlap_coefficient,&
            eligible = eligible_f,&
            eligible_by_stop_condition = eligible_by_stop_condition_f,&
            eligible_by_dimension = eligible_by_dimension_f,&
            eligible_by_var_explained = eligible_by_var_explained_f,&
            ierr = ierr&
        )

        eligible = eligible_f
        eligible_by_stop_condition = eligible_by_stop_condition_f
        eligible_by_dimension = eligible_by_dimension_f
        eligible_by_var_explained = eligible_by_var_explained_f
    end subroutine ensemble_reconciliation_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_reconciliation(module):merge_to_super_ensembles(subroutine)]]
    !| Detecting an intersection at all already requires $|\mathcal{E}_i \cap \mathcal{E}_j|$,
    !| needed by every mode; the Overlap Coefficient itself is a single extra $O(1)$ step per
    !| pair once each ensemble's own size is known --
    !| $\text{OC}(\mathcal{E}_i, \mathcal{E}_j) = |\mathcal{E}_i \cap \mathcal{E}_j| /
    !| \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$, cheaper even than the Jaccard Similarity Index
    !| it replaces (no union to derive, just the smaller of the two already-precomputed sizes)
    !| -- but modes 1 and 3 do not need it for their own decision, so its computation is
    !| guarded behind `report_overlap_coefficient`, never unconditional (see `misc/mod_STC.md`'s
    !| explicit note on this).
    !|
    !| An ineligible ensemble (`.not. eligible(i)`, see `filter_ensembles_kernel`) never
    !| contributes a pair here at all -- a plain `eligible(i) .and. eligible(j)` guard, no array
    !| copying/compaction, the same "logical AND of masks" shape mode
    !| [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]]'s
    !| own Overlap Coefficient threshold check already uses.
    !|
    !| Modes 2 and 3 group via a union-find over the qualifying-edge graph (`stc_uf_find`/
    !| `stc_uf_union` below), unioning the smaller index under the larger's root so that a
    !| component's root is always its own smallest member -- which is what makes the single
    !| pass `do r = 1, n_ensembles` below both find every component exactly once and emit them
    !| in ascending order of each group's smallest member, with no separate bookkeeping. A
    !| discovered component larger than `max_group_size` is a genuine runtime condition no
    !| static input check could foresee (it depends on the actual intersection pattern), so it
    !| is reported via `ierr` rather than silently truncated -- see `codegen_guide.md` section
    !| 5.14.
    subroutine merge_to_super_ensembles_c(&
            ensemble_masks,&
            eligible,&
            n_vectors,&
            n_ensembles,&
            mode,&
            min_overlap_coefficient,&
            report_overlap_coefficient,&
            max_group_size,&
            super_ensembles,&
            n_super_ensembles,&
            super_ensembles_overlap_coefficient,&
            ierr&
        ) bind(C, name="merge_to_super_ensembles_c")
        use tox_shape_truthful_clustering_reconciliation, only: merge_to_super_ensembles
        use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
            !! can hold; sizes its row dimension. `misc/mod_STC.md` suggests
            !! $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
            !! optional with an auto-applied default here, for the same reason as
            !! `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
            !! possibly-absent optional dummy, and a runtime-dependent value like
            !! $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
            !! default would need to be either.
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_ensembles`.
        logical(c_bool), dimension(n_vectors, n_ensembles), intent(in), target :: ensemble_masks
            !! Per-ensemble membership, see Ensemble Identification's merged output
        logical(c_bool), dimension(n_ensembles), intent(in), target :: eligible
            !! Per-ensemble eligibility to contribute a pair here at all -- see
            !! `tox_shape_truthful_clustering_filter_kernel`'s own `filter_ensembles_kernel`,
            !! this kernel's own sibling in `ensemble_reconciliation`'s two-call orchestration
        character(len=1, kind=c_char), dimension(25), intent(in), target :: mode
            !! How intersections are processed
            !!
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
            !! The default value is `1_int32`.
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical(c_bool), intent(in), target :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see the note above on this being guarded, not unconditional
            !! The default value is `.false.`.
        integer(c_int), dimension(max_group_size, n_ensembles*(n_ensembles-1)), intent(out), target :: super_ensembles
            !! One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            !! belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            !! the group's actual size, and 0 in every row of an unused trailing column beyond
            !! `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]'s
            !! own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            !! intersects) -- deliberately not divided by 2, see `ensemble_reconciliation_kernel`'s
            !! own identical note.
        integer(c_int), intent(out), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            !! actually filled
        real(c_double), dimension(max_group_size-1, n_ensembles*(n_ensembles-1)), intent(out), target :: super_ensembles_overlap_coefficient
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.
        logical, dimension(n_vectors, n_ensembles) :: ensemble_masks_f
        logical, dimension(n_ensembles) :: eligible_f
        integer(int32) :: mode_mode_f
        logical :: report_overlap_coefficient_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_NON_NULL(report_overlap_coefficient)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(mode, 25)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_ensembles*(n_ensembles-1)))
        M_CHECK_ARRAY_NON_NULL(super_ensembles_overlap_coefficient, (max_group_size-1) * (n_ensembles*(n_ensembles-1)))

        ensemble_masks_f = ensemble_masks
        eligible_f = eligible
        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("report")
                    mode_mode_f = MODE_REPORT
                case ("merge_overlap_coefficient")
                    mode_mode_f = MODE_MERGE_OVERLAP_COEFFICIENT
                case ("merge_any")
                    mode_mode_f = MODE_MERGE_ANY
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        report_overlap_coefficient_f = report_overlap_coefficient

        call merge_to_super_ensembles(&
            ensemble_masks = ensemble_masks_f,&
            eligible = eligible_f,&
            n_vectors = n_vectors,&
            n_ensembles = n_ensembles,&
            mode = mode_mode_f,&
            min_overlap_coefficient = min_overlap_coefficient,&
            report_overlap_coefficient = report_overlap_coefficient_f,&
            max_group_size = max_group_size,&
            super_ensembles = super_ensembles,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles_overlap_coefficient = super_ensembles_overlap_coefficient,&
            ierr = ierr&
        )
    end subroutine merge_to_super_ensembles_c

end module tox_shape_truthful_clustering_reconciliation_c
#endif
