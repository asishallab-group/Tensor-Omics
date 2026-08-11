#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_stc_csv(module)]]
!| Plain-text (CSV/TSV) companions to `tox_stc_json`'s JSON/HTML report, for the data-science
!| workflow the JSON alone does not serve well: loading STC's results back into a Python/R
!| analysis pipeline (`pandas.read_csv`/`read.csv`) without a JSON parser, joined back to the
!| caller's own input table purely by row number. Three independent, narrow writers, one per
!| artifact -- deliberately not derived from `tox_stc_json`'s own internal tree-building (that
!| machinery is JSON-object-shaped and private to one subroutine call, see
!| `tox_stc_json::stc_build_and_serialize_json`'s own doc comment); each writer here instead
!| recomputes its own small amount of derived state directly from the raw `ensemble_masks`/
!| `super_ensembles` arrays, the same source-of-truth `tox_stc_json` reads.
!|
!| - `serialize_stc_points_as_csv`: one row per input vector, ensemble/super-ensemble/
!| low-confidence/seed-of membership as quoted, comma-joined list cells (standard CSV
!| quoting, so `pandas.read_csv`/`read.csv` parse the embedded commas correctly) -- the
!| "assign each data point back to its ensemble(s) using the input table's own row numbers"
!| table.
!| - `serialize_stc_ensemble_overlap_as_csv`: the full pairwise Overlap Coefficient matrix
!| (only pairs with a nonempty intersection, matching `tox_stc_json`'s own convention) as a
!| plain three-column CSV.
!| - `serialize_stc_super_ensembles_as_tsv`: one line per super-ensemble, in the same
!| `<id>` TAB `<comma-separated member list>` shape as a gene-family file -- deliberately
!| *not* CSV-quoted, since the field separator (TAB) and the list separator (comma) never
!| collide.
module tox_stc_csv_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: serialize_stc_points_as_csv_c
    public :: serialize_stc_ensemble_overlap_as_csv_c
    public :: serialize_stc_super_ensembles_as_tsv_c

contains

    !> summary: C-wrapper for [[tox_stc_csv(module):serialize_stc_points_as_csv(subroutine)]]
    !| `row` is the vector's own 1-based position in the input table, so this file joins back
    !| to the caller's original data purely by row number. See this module's own doc comment
    !| for the quoting convention.
    subroutine serialize_stc_points_as_csv_c(&
            filename,&
            filename_strlen,&
            n_vectors,&
            n_selected_seed,&
            max_group_size,&
            n_super_ensembles,&
            seed_selection_mask,&
            ensemble_masks,&
            ensemble_low_confidence_masks,&
            super_ensembles,&
            ierr&
        ) bind(C, name="serialize_stc_points_as_csv_c")
        use tox_stc_csv, only: serialize_stc_points_as_csv

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the CSV file to write
        integer(c_int), intent(in), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        logical(c_bool), dimension(n_vectors), intent(in), target :: seed_selection_mask
            !! Seed selection, see `seeds`
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership
        integer(c_int), dimension(max_group_size, n_selected_seed*(n_selected_seed-1)), intent(in), target :: super_ensembles
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        character(len=:), allocatable :: filename_f
        logical, dimension(n_vectors) :: seed_selection_mask_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_masks_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_low_confidence_masks_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_seed)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(seed_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_selected_seed*(n_selected_seed-1)))

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        seed_selection_mask_f = seed_selection_mask
        ensemble_masks_f = ensemble_masks
        ensemble_low_confidence_masks_f = ensemble_low_confidence_masks

        call serialize_stc_points_as_csv(&
            filename = filename_f,&
            n_vectors = n_vectors,&
            n_selected_seed = n_selected_seed,&
            max_group_size = max_group_size,&
            n_super_ensembles = n_super_ensembles,&
            seed_selection_mask = seed_selection_mask_f,&
            ensemble_masks = ensemble_masks_f,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks_f,&
            super_ensembles = super_ensembles,&
            ierr = ierr&
        )
    end subroutine serialize_stc_points_as_csv_c

    !> summary: C-wrapper for [[tox_stc_csv(module):serialize_stc_ensemble_overlap_as_csv(subroutine)]]
    !| Only pairs with a nonempty intersection are written, matching `tox_stc_json`'s own
    !| `overlap_coefficient_matrix` convention -- an absent pair means Overlap Coefficient 0.
    subroutine serialize_stc_ensemble_overlap_as_csv_c(&
            filename,&
            filename_strlen,&
            n_vectors,&
            n_selected_seed,&
            ensemble_masks,&
            ierr&
        ) bind(C, name="serialize_stc_ensemble_overlap_as_csv_c")
        use tox_stc_csv, only: serialize_stc_ensemble_overlap_as_csv

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the CSV file to write
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        character(len=:), allocatable :: filename_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_masks_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_selected_seed)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        ensemble_masks_f = ensemble_masks

        call serialize_stc_ensemble_overlap_as_csv(&
            filename = filename_f,&
            n_vectors = n_vectors,&
            n_selected_seed = n_selected_seed,&
            ensemble_masks = ensemble_masks_f,&
            ierr = ierr&
        )
    end subroutine serialize_stc_ensemble_overlap_as_csv_c

    !> summary: C-wrapper for [[tox_stc_csv(module):serialize_stc_super_ensembles_as_tsv(subroutine)]]
    !| One line per super-ensemble: `<group_id>` TAB `<comma-separated member ensemble ids>`,
    !| no header, no quoting -- see this module's own doc comment for why.
    subroutine serialize_stc_super_ensembles_as_tsv_c(&
            filename,&
            filename_strlen,&
            max_group_size,&
            n_super_ensembles,&
            super_ensembles,&
            ierr&
        ) bind(C, name="serialize_stc_super_ensembles_as_tsv_c")
        use tox_stc_csv, only: serialize_stc_super_ensembles_as_tsv

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(c_int), intent(in), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the TSV file to write
        integer(c_int), dimension(max_group_size, n_super_ensembles), intent(in), target :: super_ensembles
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        character(len=:), allocatable :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * n_super_ensembles)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return

        call serialize_stc_super_ensembles_as_tsv(&
            filename = filename_f,&
            max_group_size = max_group_size,&
            n_super_ensembles = n_super_ensembles,&
            super_ensembles = super_ensembles,&
            ierr = ierr&
        )
    end subroutine serialize_stc_super_ensembles_as_tsv_c

end module tox_stc_csv_c
#endif
