#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Ensemble Reconciliation
!|
!| `ensemble_reconciliation`: identifies intersecting ensembles from Ensemble Identification's
!| merged `ensemble_masks` output and, depending on `mode`, either just reports intersecting
!| pairs or groups transitively-intersecting ensembles into "super-ensembles" via a union-find
!| over the pairwise intersection graph. See `misc/mod_STC.md`, "Ensemble Reconciliation", for
!| the full algorithm definition. Does not alter Ensemble Identification's own result -- this
!| module only reports and groups, on the side.
!|
!| A thin, two-call orchestrator over its own two sibling kernels: first
!| [[tox_shape_truthful_clustering_filter_impl(module):filter_ensembles_impl]] (which
!| ensembles are even eligible to contribute a pair, by Stop Condition/final dimension/final
!| variance explained), then this module's own `merge_to_super_ensembles_impl` (the actual
!| pairwise-intersection/union-find grouping, over eligible ensembles only). Splitting these
!| into two independently testable, independently reusable kernels -- rather than one kernel
!| that both decides eligibility and merges -- is a deliberate design choice: eligibility is a
!| statement about *individual* ensembles (their own Stop Condition/geometry), merging is a
!| statement about *pairs*, and conflating the two made every new filtering criterion require
!| touching the same monolithic merge logic. See `misc/mod_STC.md`'s own rationale for the
!| split.
module tox_shape_truthful_clustering_reconciliation_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, ERR_SIZE_MISMATCH
    use tox_shape_truthful_clustering_filter_impl, only: filter_ensembles_impl
    M_IMPLICIT_NONE

#define CM_STC_MODE_REPORT 1_int32
#define CM_STC_MODE_MERGE_OVERLAP_COEFFICIENT 2_int32
#define CM_STC_MODE_MERGE_ANY 3_int32
#define CM_STC_MIN_OVERLAP_COEFFICIENT_DEFAULT 0.9_real64

    private

    !> Mode 1: report every intersecting pair as its own 2-row column, no transitive grouping.
    integer(int32), parameter, public :: MODE_REPORT = CM_STC_MODE_REPORT
    !> Mode 2: transitively group ensembles connected by an edge whose Overlap Coefficient is
    !| >= `min_overlap_coefficient`.
    integer(int32), parameter, public :: MODE_MERGE_OVERLAP_COEFFICIENT = CM_STC_MODE_MERGE_OVERLAP_COEFFICIENT
    !> Mode 3: transitively group ensembles connected by any nonempty intersection.
    integer(int32), parameter, public :: MODE_MERGE_ANY = CM_STC_MODE_MERGE_ANY

    public :: ensemble_reconciliation_impl
    public :: merge_to_super_ensembles_impl

contains

    !> summary: Filter eligible ensembles, then group/report their intersections
    !| AUTHOR_ASIS_HALLAB
    !| See this module's own header comment for why this is a two-call orchestrator, and
    !| `tox_shape_truthful_clustering_filter_impl`'s own `filter_ensembles_impl` for the
    !| full eligibility-filtering algorithm (Stop Condition/final dimension/final variance
    !| explained, each independently optional, combined by logical AND). `ierr` is set only if
    !| `merge_to_super_ensembles_impl` discovers a component larger than `max_group_size` --
    !| see that kernel's own doc comment.
    pure subroutine ensemble_reconciliation_impl(ensemble_masks, ensemble_stop_reason, n_dimensions, n_vectors, &
                                                    n_ensembles, &
                                                    ensemble_U_history, ensemble_d_history, ensemble_S_history, &
                                                    ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                                                    ensemble_accepted_history, o, &
                                                    mode, &
                                                    min_overlap_coefficient, report_overlap_coefficient, &
                                                    allowed_stop_reasons, filter_dim_min, filter_dim_max, var_explained_min, &
                                                    max_group_size, &
                                                    super_ensembles, n_super_ensembles, &
                                                    super_ensembles_overlap_coefficient, &
                                                    eligible, eligible_by_stop_condition, eligible_by_dimension, &
                                                    eligible_by_var_explained, &
                                                    ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E, see Ensemble Identification's merged `ensemble_masks`.
            !! At least 2: with fewer, no pair can ever intersect, so there is nothing this
            !! module could report or group.
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth
            !! DM_MIN(1_int32)
        logical(c_bool), intent(in) :: ensemble_masks(n_vectors, n_ensembles)
            !! Per-ensemble membership, see Ensemble Identification's merged output
        integer(int32), intent(in) :: ensemble_stop_reason(n_ensembles)
            !! Per-ensemble Stop Condition, see `filter_ensembles_impl`
            !! DM_MIN(1_int32)
            !! DM_MAX(4_int32)
        real(real64), intent(in) :: ensemble_U_history(n_dimensions, n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing tangent+normal bases, see Ensemble Identification's merged
            !! output
        integer(int32), intent(in) :: ensemble_d_history(o, n_ensembles)
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), intent(in) :: ensemble_S_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing singular values
        real(real64), intent(in) :: ensemble_mu_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing centers
        real(real64), intent(in) :: ensemble_G_history(o, n_ensembles)
            !! Per-ensemble trailing spectral gaps
        integer(int32), intent(in) :: ensemble_k_history(o, n_ensembles)
            !! Per-ensemble trailing sizes
        logical(c_bool), intent(in) :: ensemble_accepted_history(o, n_ensembles)
            !! Whether the growth iteration retained in each history column was itself accepted
        integer(int32), intent(in), optional :: mode
            !! How intersections are processed
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Report intersecting pairs only        | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]     |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]]  |
            !! | Merge transitively on any intersection | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_ANY(variable)]]  |
            !! DM_DEFAULT(CM_STC_MODE_REPORT)
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_STC_MIN_OVERLAP_COEFFICIENT_DEFAULT)
        logical(c_bool), intent(in), optional :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see `merge_to_super_ensembles_impl`'s own note on this being guarded, not
            !! unconditional
            !! DM_DEFAULT(.false.)
        logical(c_bool), intent(in), optional :: allowed_stop_reasons(4)
            !! See `tox_shape_truthful_clustering_filter_impl`'s own
            !! `filter_ensembles_by_stop_condition_impl`
        integer(int32), intent(in), optional :: filter_dim_min
            !! See `tox_shape_truthful_clustering_filter_impl`'s own
            !! `filter_ensembles_by_dimension_impl`
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in), optional :: filter_dim_max
            !! See `tox_shape_truthful_clustering_filter_impl`'s own
            !! `filter_ensembles_by_dimension_impl`
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in), optional :: var_explained_min
            !! See `tox_shape_truthful_clustering_filter_impl`'s own
            !! `filter_ensembles_by_var_explained_impl`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), intent(in) :: max_group_size
            !! Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
            !! can hold; sizes its row dimension. `misc/mod_STC.md` suggests
            !! $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
            !! optional with an auto-applied default here, for the same reason as
            !! `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
            !! possibly-absent optional dummy, and a runtime-dependent value like
            !! $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
            !! default would need to be either.
            !! DM_MIN(2_int32)
            !! DM_MAX(n_ensembles)
        integer(int32), intent(out) :: super_ensembles(max_group_size, n_ensembles*(n_ensembles - 1))
            !! One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            !! belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            !! the group's actual size, and 0 in every row of an unused trailing column beyond
            !! `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            !! [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]'s
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
        real(real64), intent(out) :: super_ensembles_overlap_coefficient(max_group_size - 1, n_ensembles*(n_ensembles - 1))
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        logical(c_bool), intent(out) :: eligible(n_ensembles)
            !! Combined per-ensemble eligibility actually used for merging above -- see
            !! `filter_ensembles_impl`. Ineligible ensembles are otherwise untouched: they
            !! are never removed from `ensemble_masks` or anything else this whole family
            !! reports, only excluded from contributing a pair here.
        logical(c_bool), intent(out) :: eligible_by_stop_condition(n_ensembles)
            !! See `filter_ensembles_by_stop_condition_impl`
        logical(c_bool), intent(out) :: eligible_by_dimension(n_ensembles)
            !! See `filter_ensembles_by_dimension_impl`
        logical(c_bool), intent(out) :: eligible_by_var_explained(n_ensembles)
            !! See `filter_ensembles_by_var_explained_impl`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.

        call set_ok(ierr)

        call filter_ensembles_impl(n_dimensions, o, n_ensembles, &
            ensemble_U_history, ensemble_d_history, ensemble_S_history, &
            ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, &
            ensemble_stop_reason, allowed_stop_reasons, filter_dim_min, filter_dim_max, var_explained_min, &
            eligible, eligible_by_stop_condition, eligible_by_dimension, eligible_by_var_explained)

        call merge_to_super_ensembles_impl(ensemble_masks, eligible, n_vectors, n_ensembles, &
            mode, min_overlap_coefficient, report_overlap_coefficient, max_group_size, &
            super_ensembles, n_super_ensembles, super_ensembles_overlap_coefficient, ierr)

    end subroutine ensemble_reconciliation_impl

    !> summary: Group/report intersections among eligible ensembles into super-ensembles
    !| AUTHOR_ASIS_HALLAB
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
    !| An ineligible ensemble (`.not. eligible(i)`, see `filter_ensembles_impl`) never
    !| contributes a pair here at all -- a plain `eligible(i) .and. eligible(j)` guard, no array
    !| copying/compaction, the same "logical AND of masks" shape mode
    !| [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]]'s
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
    pure subroutine merge_to_super_ensembles_impl(ensemble_masks, eligible, n_vectors, n_ensembles, &
                                                     mode, min_overlap_coefficient, report_overlap_coefficient, &
                                                     max_group_size, &
                                                     super_ensembles, n_super_ensembles, &
                                                     super_ensembles_overlap_coefficient, ierr)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! DM_MIN(2_int32)
        logical(c_bool), intent(in) :: ensemble_masks(n_vectors, n_ensembles)
            !! Per-ensemble membership, see Ensemble Identification's merged output
        logical(c_bool), intent(in) :: eligible(n_ensembles)
            !! Per-ensemble eligibility to contribute a pair here at all -- see
            !! `tox_shape_truthful_clustering_filter_impl`'s own `filter_ensembles_impl`,
            !! this kernel's own sibling in `ensemble_reconciliation`'s two-call orchestration
        integer(int32), intent(in), optional :: mode
            !! How intersections are processed
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Report intersecting pairs only        | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]     |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]]  |
            !! | Merge transitively on any intersection | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_ANY(variable)]]  |
            !! DM_DEFAULT(CM_STC_MODE_REPORT)
        real(real64), intent(in), optional :: min_overlap_coefficient
            !! Minimum Overlap Coefficient ($|\mathcal{E}_i \cap \mathcal{E}_j| /
            !! \min(|\mathcal{E}_i|, |\mathcal{E}_j|)$) for an edge to qualify in mode
            !! [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]];
            !! ignored in every other mode
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_STC_MIN_OVERLAP_COEFFICIENT_DEFAULT)
        logical(c_bool), intent(in), optional :: report_overlap_coefficient
            !! Whether to compute and return `super_ensembles_overlap_coefficient` at all --
            !! see the note above on this being guarded, not unconditional
            !! DM_DEFAULT(.false.)
        integer(int32), intent(in) :: max_group_size
            !! Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
            !! can hold; sizes its row dimension. `misc/mod_STC.md` suggests
            !! $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
            !! optional with an auto-applied default here, for the same reason as
            !! `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
            !! possibly-absent optional dummy, and a runtime-dependent value like
            !! $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
            !! default would need to be either.
            !! DM_MIN(2_int32)
            !! DM_MAX(n_ensembles)
        integer(int32), intent(out) :: super_ensembles(max_group_size, n_ensembles*(n_ensembles - 1))
            !! One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
            !! belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
            !! the group's actual size, and 0 in every row of an unused trailing column beyond
            !! `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
            !! [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]'s
            !! own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
            !! intersects) -- deliberately not divided by 2, see `ensemble_reconciliation_impl`'s
            !! own identical note.
        integer(int32), intent(out) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles`/`super_ensembles_overlap_coefficient`
            !! actually filled
        real(real64), intent(out) :: super_ensembles_overlap_coefficient(max_group_size - 1, n_ensembles*(n_ensembles - 1))
            !! Column $l$, row $c_i$: the Overlap Coefficient between the ensembles in
            !! `super_ensembles(c_i, l)` and `super_ensembles(c_i + 1, l)`. All zero unless
            !! `report_overlap_coefficient` was requested -- see the note above.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success. Set only if a discovered component's size exceeds
            !! `max_group_size` -- not a condition any input check could foresee, see above.

        integer(int32) :: ensemble_size(n_ensembles)
        integer(int32) :: parent(n_ensembles)
        integer(int32) :: member(n_ensembles)
        integer(int32) :: actual_mode
        real(real64)   :: actual_min_overlap_coefficient
        logical(c_bool)        :: actual_report_overlap_coefficient
        integer(int32) :: i, j, l, r, group_size, intersect_count, root

        call set_ok(ierr)

        M_DEFAULT_VAL(mode, actual_mode, CM_STC_MODE_REPORT)
        M_DEFAULT_VAL(min_overlap_coefficient, actual_min_overlap_coefficient, CM_STC_MIN_OVERLAP_COEFFICIENT_DEFAULT)
        M_DEFAULT_VAL(report_overlap_coefficient, actual_report_overlap_coefficient, .false.)

        super_ensembles                       = 0
        super_ensembles_overlap_coefficient   = 0.0_real64
        n_super_ensembles                     = 0

        do i = 1, n_ensembles
            ensemble_size(i) = count(ensemble_masks(:, i))
            parent(i)        = i
        end do

        if (actual_mode == CM_STC_MODE_REPORT) then

            do i = 1, n_ensembles - 1
                do j = i + 1, n_ensembles
                    if (.not. eligible(i) .or. .not. eligible(j)) cycle
                    intersect_count = count(ensemble_masks(:, i) .and. ensemble_masks(:, j))
                    if (intersect_count < 1) cycle

                    n_super_ensembles = n_super_ensembles + 1
                    super_ensembles(1, n_super_ensembles) = i
                    super_ensembles(2, n_super_ensembles) = j
                    if (actual_report_overlap_coefficient) then
                        super_ensembles_overlap_coefficient(1, n_super_ensembles) = real(intersect_count, real64) / &
                            real(min(ensemble_size(i), ensemble_size(j)), real64)
                    end if
                end do
            end do

        else

            do i = 1, n_ensembles - 1
                do j = i + 1, n_ensembles
                    if (.not. eligible(i) .or. .not. eligible(j)) cycle
                    intersect_count = count(ensemble_masks(:, i) .and. ensemble_masks(:, j))
                    if (intersect_count < 1) cycle
                    if (actual_mode == CM_STC_MODE_MERGE_OVERLAP_COEFFICIENT) then
                        if (real(intersect_count, real64) / &
                            real(min(ensemble_size(i), ensemble_size(j)), real64) < actual_min_overlap_coefficient) cycle
                    end if
                    call stc_uf_union(parent, n_ensembles, i, j)
                end do
            end do

            do i = 1, n_ensembles
                call stc_uf_find(parent, n_ensembles, i, root)
                parent(i) = root
            end do

            do r = 1, n_ensembles
                group_size = count(parent == r)
                if (group_size < 2) cycle

                if (group_size > max_group_size) then
                    call set_err_once(ierr, ERR_SIZE_MISMATCH)
                    return
                end if

                l = 0
                do i = 1, n_ensembles
                    if (parent(i) /= r) cycle
                    l = l + 1
                    member(l) = i
                end do

                n_super_ensembles = n_super_ensembles + 1
                super_ensembles(1:group_size, n_super_ensembles) = member(1:group_size)

                if (actual_report_overlap_coefficient) then
                    do l = 1, group_size - 1
                        intersect_count = count(ensemble_masks(:, member(l)) .and. ensemble_masks(:, member(l + 1)))
                        super_ensembles_overlap_coefficient(l, n_super_ensembles) = real(intersect_count, real64) / &
                            real(min(ensemble_size(member(l)), ensemble_size(member(l + 1))), real64)
                    end do
                end if
            end do

        end if

    end subroutine merge_to_super_ensembles_impl

    !> Union-find root lookup with path compression. A `pure` function's dummy arguments
    !| must all be `intent(in)` -- path compression mutates `parent`, so this has to be a
    !| subroutine, not a function. Not itself a kernel (private, no wrapper): shared
    !| bookkeeping for `merge_to_super_ensembles_impl`, not part of the public contract.
    recursive pure subroutine stc_uf_find(parent, n, i, root)
        integer(int32), intent(in) :: n
        integer(int32), intent(inout) :: parent(n)
        integer(int32), intent(in) :: i
        integer(int32), intent(out) :: root

        if (parent(i) == i) then
            root = i
        else
            call stc_uf_find(parent, n, parent(i), root)
            parent(i) = root
        end if
    end subroutine stc_uf_find

    !> Union-find merge, always attaching the numerically larger root under the smaller --
    !| so a component's root is always its own smallest member, see
    !| `merge_to_super_ensembles_impl`'s own doc comment above.
    pure subroutine stc_uf_union(parent, n, i, j)
        integer(int32), intent(in) :: n
        integer(int32), intent(inout) :: parent(n)
        integer(int32), intent(in) :: i, j
        integer(int32) :: root_i, root_j

        call stc_uf_find(parent, n, i, root_i)
        call stc_uf_find(parent, n, j, root_j)
        if (root_i == root_j) return

        if (root_i < root_j) then
            parent(root_j) = root_i
        else
            parent(root_i) = root_j
        end if
    end subroutine stc_uf_union

end module tox_shape_truthful_clustering_reconciliation_impl
