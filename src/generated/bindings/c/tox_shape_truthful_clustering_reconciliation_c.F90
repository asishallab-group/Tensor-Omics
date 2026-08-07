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

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_reconciliation(module):ensemble_reconciliation(subroutine)]]
    !| Detecting an intersection at all already requires $|\mathcal{E}_i \cap \mathcal{E}_j|$,
    !| needed by every mode; the JSI itself is a single extra $O(1)$ step per pair once each
    !| ensemble's own size is known, via
    !| $|\mathcal{E}_i \cup \mathcal{E}_j| = |\mathcal{E}_i| + |\mathcal{E}_j| - |\mathcal{E}_i \cap \mathcal{E}_j|$
    !| -- but modes 1 and 3 do not need it for their own decision, so its computation is
    !| guarded behind `report_jsi`, never unconditional (see `misc/mod_STC.md`'s explicit note
    !| on this).
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
    subroutine ensemble_reconciliation_c(&
            ensemble_masks,&
            n_vectors,&
            n_ensembles,&
            mode,&
            min_jsi,&
            report_jsi,&
            max_group_size,&
            super_ensembles,&
            n_super_ensembles,&
            super_ensembles_JSI,&
            ierr&
        ) bind(C, name="ensemble_reconciliation_c")
        use tox_shape_truthful_clustering_reconciliation, only: ensemble_reconciliation
        use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_JSI, MODE_REPORT

        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E, see Ensemble Identification's merged `ensemble_masks`.
            !! At least 2: with fewer, no pair can ever intersect, so there is nothing this
            !! module could report or group.
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
        character(len=1, kind=c_char), dimension(9), intent(in), target :: mode
            !! How intersections are processed
            !!
            !! | Mode                                   | Value                                                                                    |
            !! |----------------------------------------|------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only         | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]    |
            !! | Merge transitively at a minimum JSI    | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_JSI(variable)]] |
            !! | Merge transitively on any intersection | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]] |
            !! The default value is `1_int32`.
        real(c_double), intent(in), target :: min_jsi
            !! Minimum Jaccard Similarity Index for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_JSI(variable)]];
            !! ignored in every other mode
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        logical(c_bool), intent(in), target :: report_jsi
            !! Whether to compute and return `super_ensembles_JSI` at all -- see the note
            !! above on this being guarded, not unconditional
            !! The default value is `.false.`.
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
            !! Number of leading columns of `super_ensembles`/`super_ensembles_JSI` actually
            !! filled
        real(c_double), dimension(max_group_size-1, n_ensembles*(n_ensembles-1)), intent(out), target :: super_ensembles_JSI
            !! Column $l$, row $c_i$: the JSI between the ensembles in `super_ensembles(c_i, l)`
            !! and `super_ensembles(c_i + 1, l)`. All zero unless `report_jsi` was requested --
            !! see the note above.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.
        logical, dimension(n_vectors, n_ensembles) :: ensemble_masks_f
        integer(int32) :: mode_mode_f
        logical :: report_jsi_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_NON_NULL(min_jsi)
        M_CHECK_NON_NULL(report_jsi)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(mode, 9)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_ensembles*(n_ensembles-1)))
        M_CHECK_ARRAY_NON_NULL(super_ensembles_JSI, (max_group_size-1) * (n_ensembles*(n_ensembles-1)))

        ensemble_masks_f = ensemble_masks
        block
            character(len=:), allocatable :: mode_f
            call c_char_1d_as_string(mode, mode_f, ierr)
            if (is_err(ierr)) return

            select case (mode_f)
                case ("report")
                    mode_mode_f = MODE_REPORT
                case ("merge_jsi")
                    mode_mode_f = MODE_MERGE_JSI
                case ("merge_any")
                    mode_mode_f = MODE_MERGE_ANY
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        report_jsi_f = report_jsi

        call ensemble_reconciliation(&
            ensemble_masks = ensemble_masks_f,&
            n_vectors = n_vectors,&
            n_ensembles = n_ensembles,&
            mode = mode_mode_f,&
            min_jsi = min_jsi,&
            report_jsi = report_jsi_f,&
            max_group_size = max_group_size,&
            super_ensembles = super_ensembles,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles_JSI = super_ensembles_JSI,&
            ierr = ierr&
        )
    end subroutine ensemble_reconciliation_c

end module tox_shape_truthful_clustering_reconciliation_c
#endif
