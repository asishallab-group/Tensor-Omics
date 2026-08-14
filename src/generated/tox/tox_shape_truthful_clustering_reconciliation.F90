#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_reconciliation_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_reconciliation
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT, ensemble_reconciliation_kernel
    use tox_shape_truthful_clustering_reconciliation_kernel, only: merge_to_super_ensembles_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, clear_err_arg_pos
    use tox_errors, only: set_err_once, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: ensemble_reconciliation
    public :: merge_to_super_ensembles

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_reconciliation_kernel(module):ensemble_reconciliation_kernel]].
    !| See this module's own header comment for why this is a two-call orchestrator, and
    !| `tox_shape_truthful_clustering_filter_kernel`'s own `filter_ensembles_kernel` for the
    !| full eligibility-filtering algorithm (Stop Condition/final dimension/final variance
    !| explained, each independently optional, combined by logical AND). `ierr` is set only if
    !| `merge_to_super_ensembles_kernel` discovers a component larger than `max_group_size` --
    !| see that kernel's own doc comment.
    subroutine ensemble_reconciliation(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E, see Ensemble Identification's merged `ensemble_masks`.
            !! At least 2: with fewer, no pair can ever intersect, so there is nothing this
            !! module could report or group.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_group_size
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
        logical, dimension(n_vectors, n_ensembles), intent(in) :: ensemble_masks
            !! Per-ensemble membership, see Ensemble Identification's merged output
        integer(int32), dimension(n_ensembles), intent(in) :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `filter_ensembles_kernel`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `4_int32`.
        real(real64), dimension(n_dimensions, n_dimensions, o, n_ensembles), intent(in) :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see Ensemble Identification's merged
            !! output
        integer(int32), dimension(o, n_ensembles), intent(in) :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), dimension(n_dimensions, o, n_ensembles), intent(in) :: ensemble_S_history
            !! Per-ensemble trailing singular values
        real(real64), dimension(n_dimensions, o, n_ensembles), intent(in) :: ensemble_mu_history
            !! Per-ensemble trailing centers
        real(real64), dimension(o, n_ensembles), intent(in) :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        integer(int32), dimension(o, n_ensembles), intent(in) :: ensemble_k_history
            !! Per-ensemble trailing sizes
        logical, dimension(o, n_ensembles), intent(in) :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself accepted
        integer(int32), intent(in), optional :: mode
            !! How intersections are processed
            !!
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical, intent(in), optional :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see `merge_to_super_ensembles_kernel`'s own note on this being guarded, not
            !! unconditional
            !! The default value is `.false.`.
        logical, dimension(4), intent(in), optional :: allowed_stop_reasons
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_stop_condition_kernel`
        integer(int32), intent(in), optional :: d_min
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_dimension_kernel`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in), optional :: d_max
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_dimension_kernel`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), intent(in), optional :: var_explained_min
            !! See `tox_shape_truthful_clustering_filter_kernel`'s own
            !! `filter_ensembles_by_var_explained_kernel`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(max_group_size, n_ensembles*(n_ensembles-1)), intent(out) :: super_ensembles
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
        integer(int32), intent(out) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            !! actually filled
        real(real64), dimension(max_group_size-1, n_ensembles*(n_ensembles-1)), intent(out) :: super_ensembles_overlap_coefficient
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        logical, dimension(n_ensembles), intent(out) :: eligible
            !! Combined per-ensemble eligibility actually used for merging above -- see
            !! `filter_ensembles_kernel`. Ineligible ensembles are otherwise untouched: they
            !! are never removed from `ensemble_masks` or anything else this whole family
            !! reports, only excluded from contributing a pair here.
        logical, dimension(n_ensembles), intent(out) :: eligible_by_stop_condition
            !! See `filter_ensembles_by_stop_condition_kernel`
        logical, dimension(n_ensembles), intent(out) :: eligible_by_dimension
            !! See `filter_ensembles_by_dimension_kernel`
        logical, dimension(n_ensembles), intent(out) :: eligible_by_var_explained
            !! See `filter_ensembles_by_var_explained_kernel`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=3_int32, min=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_ensembles, ierr, arg_pos=5_int32, min=2_int32)
        call validate_in_range_int(o, ierr, arg_pos=13_int32, min=1_int32)
        call validate_in_range_real(min_overlap_coefficient, ierr, arg_pos=15_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(d_min, ierr, arg_pos=18_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_int(d_max, ierr, arg_pos=19_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(var_explained_min, ierr, arg_pos=20_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(max_group_size, ierr, arg_pos=21_int32, min=2_int32, max=n_ensembles)
        call validate_all_in_range_int(ensemble_stop_reason, n_ensembles, ierr, arg_pos=2_int32, min=1_int32, max=4_int32)
        call validate_all_in_range_real(ensemble_U_history, n_dimensions * n_dimensions * o * n_ensembles, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(ensemble_S_history, n_dimensions * o * n_ensembles, ierr, arg_pos=8_int32)
        call validate_all_in_range_real(ensemble_mu_history, n_dimensions * o * n_ensembles, ierr, arg_pos=9_int32)
        call validate_all_in_range_real(ensemble_G_history, o * n_ensembles, ierr, arg_pos=10_int32)
        if (present(mode)) then; if (mode /= MODE_REPORT .and. mode /= MODE_MERGE_OVERLAP_COEFFICIENT .and. mode /= MODE_MERGE_ANY) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=14_int32); end if
        if (is_err(ierr)) return
#endif

        call ensemble_reconciliation_kernel(&
            ensemble_masks = ensemble_masks,&
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
            ensemble_accepted_history = ensemble_accepted_history,&
            o = o,&
            mode = mode,&
            min_overlap_coefficient = min_overlap_coefficient,&
            report_overlap_coefficient = report_overlap_coefficient,&
            allowed_stop_reasons = allowed_stop_reasons,&
            d_min = d_min,&
            d_max = d_max,&
            var_explained_min = var_explained_min,&
            max_group_size = max_group_size,&
            super_ensembles = super_ensembles,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles_overlap_coefficient = super_ensembles_overlap_coefficient,&
            eligible = eligible,&
            eligible_by_stop_condition = eligible_by_stop_condition,&
            eligible_by_dimension = eligible_by_dimension,&
            eligible_by_var_explained = eligible_by_var_explained,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine ensemble_reconciliation

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_reconciliation_kernel(module):merge_to_super_ensembles_kernel]].
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
    subroutine merge_to_super_ensembles(&
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
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: max_group_size
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
        logical, dimension(n_vectors, n_ensembles), intent(in) :: ensemble_masks
            !! Per-ensemble membership, see Ensemble Identification's merged output
        logical, dimension(n_ensembles), intent(in) :: eligible
            !! Per-ensemble eligibility to contribute a pair here at all -- see
            !! `tox_shape_truthful_clustering_filter_kernel`'s own `filter_ensembles_kernel`,
            !! this kernel's own sibling in `ensemble_reconciliation`'s two-call orchestration
        integer(int32), intent(in), optional :: mode
            !! How intersections are processed
            !!
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
            !! The default value is `1_int32`.
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        logical, intent(in), optional :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see the note above on this being guarded, not unconditional
            !! The default value is `.false.`.
        integer(int32), dimension(max_group_size, n_ensembles*(n_ensembles-1)), intent(out) :: super_ensembles
            !! One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            !! belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            !! the group's actual size, and 0 in every row of an unused trailing column beyond
            !! `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]'s
            !! own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            !! intersects) -- deliberately not divided by 2, see `ensemble_reconciliation_kernel`'s
            !! own identical note.
        integer(int32), intent(out) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            !! actually filled
        real(real64), dimension(max_group_size-1, n_ensembles*(n_ensembles-1)), intent(out) :: super_ensembles_overlap_coefficient
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_ensembles, ierr, arg_pos=4_int32, min=2_int32)
        call validate_in_range_real(min_overlap_coefficient, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(max_group_size, ierr, arg_pos=8_int32, min=2_int32, max=n_ensembles)
        if (present(mode)) then; if (mode /= MODE_REPORT .and. mode /= MODE_MERGE_OVERLAP_COEFFICIENT .and. mode /= MODE_MERGE_ANY) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32); end if
        if (is_err(ierr)) return
#endif

        call merge_to_super_ensembles_kernel(&
            ensemble_masks = ensemble_masks,&
            eligible = eligible,&
            n_vectors = n_vectors,&
            n_ensembles = n_ensembles,&
            mode = mode,&
            min_overlap_coefficient = min_overlap_coefficient,&
            report_overlap_coefficient = report_overlap_coefficient,&
            max_group_size = max_group_size,&
            super_ensembles = super_ensembles,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles_overlap_coefficient = super_ensembles_overlap_coefficient,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine merge_to_super_ensembles

end module tox_shape_truthful_clustering_reconciliation
