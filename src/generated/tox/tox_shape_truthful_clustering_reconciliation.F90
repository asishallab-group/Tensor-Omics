#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_reconciliation_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_reconciliation
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT, ensemble_reconciliation_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, clear_err_arg_pos
    use tox_errors, only: set_err_once, validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: ensemble_reconciliation

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_reconciliation_kernel(module):ensemble_reconciliation_kernel]].
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
    !| Modes 2 and 3 group via a union-find over the qualifying-edge graph (`stc_uf_find`/
    !| `stc_uf_union` below), unioning the smaller index under the larger's root so that a
    !| component's root is always its own smallest member -- which is what makes the single
    !| pass `do r = 1, n_ensembles` below both find every component exactly once and emit them
    !| in ascending order of each group's smallest member, with no separate bookkeeping. A
    !| discovered component larger than `max_group_size` is a genuine runtime condition no
    !| static input check could foresee (it depends on the actual intersection pattern), so it
    !| is reported via `ierr` rather than silently truncated -- see `codegen_guide.md` section
    !| 5.14.
    subroutine ensemble_reconciliation(&
            ensemble_masks,&
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
            !! Number of ensembles N_E, see Ensemble Identification's merged `ensemble_masks`.
            !! At least 2: with fewer, no pair can ever intersect, so there is nothing this
            !! module could report or group.
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
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_ensembles, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_real(min_overlap_coefficient, ierr, arg_pos=5_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_int(max_group_size, ierr, arg_pos=7_int32, min=2_int32, max=n_ensembles)
        if (present(mode)) then; if (mode /= MODE_REPORT .and. mode /= MODE_MERGE_OVERLAP_COEFFICIENT .and. mode /= MODE_MERGE_ANY) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32); end if
        if (is_err(ierr)) return
#endif

        call ensemble_reconciliation_kernel(&
            ensemble_masks = ensemble_masks,&
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
    end subroutine ensemble_reconciliation

end module tox_shape_truthful_clustering_reconciliation
