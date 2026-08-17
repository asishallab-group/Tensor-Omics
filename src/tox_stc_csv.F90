#include <src/macros.h>

!> Plain-text (CSV/TSV) companions to `tox_stc_json`'s JSON/HTML report, for the data-science
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
!|   low-confidence/seed-of membership as quoted, comma-joined list cells (standard CSV
!|   quoting, so `pandas.read_csv`/`read.csv` parse the embedded commas correctly) -- the
!|   "assign each data point back to its ensemble(s) using the input table's own row numbers"
!|   table.
!| - `serialize_stc_ensemble_overlap_as_csv`: the full pairwise Overlap Coefficient matrix
!|   (only pairs with a nonempty intersection, matching `tox_stc_json`'s own convention) as a
!|   plain three-column CSV.
!| - `serialize_stc_super_ensembles_as_tsv`: one line per super-ensemble, in the same
!|   `<id>` TAB `<comma-separated member list>` shape as a gene-family file -- deliberately
!|   *not* CSV-quoted, since the field separator (TAB) and the list separator (comma) never
!|   collide.
module tox_stc_csv
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    use f42_safeguard
    use f42_io, only: open_file
    use tox_errors, only: is_err, set_ok, validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE

    private
    public :: serialize_stc_points_as_csv
    public :: serialize_stc_ensemble_overlap_as_csv
    public :: serialize_stc_super_ensembles_as_tsv

contains

    !> M_EXPORT_C
    !| summary: Serializes each input vector's ensemble/super-ensemble membership as CSV
    !| AUTHOR_ASIS_HALLAB
    !| `row` is the vector's own 1-based position in the input table, so this file joins back
    !| to the caller's original data purely by row number. See this module's own doc comment
    !| for the quoting convention.
    subroutine serialize_stc_points_as_csv(filename, n_vectors, n_selected_seed, max_group_size, &
                                           n_super_ensembles, seed_selection_mask, ensemble_masks, &
                                           ensemble_low_confidence_masks, super_ensembles, ierr)
        character(len=*), intent(in) :: filename
            !! Name of the CSV file to write
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(int32), intent(in) :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(int32), intent(in) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        logical(c_bool), intent(in) :: seed_selection_mask(n_vectors)
            !! Seed selection, see `seeds`
        logical(c_bool), intent(in) :: ensemble_masks(n_vectors, n_selected_seed)
            !! Per-ensemble accepted membership, one column per seed
        logical(c_bool), intent(in) :: ensemble_low_confidence_masks(n_vectors, n_selected_seed)
            !! Per-ensemble iteration-1 fallback membership
        integer(int32), intent(in) :: super_ensembles(max_group_size, n_selected_seed*(n_selected_seed - 1))
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: ensemble_super_id(n_selected_seed)
        integer(int32) :: vector_to_seed_index(n_vectors)
        integer(int32) :: ensembles_of_point(n_selected_seed)
        integer(int32) :: super_ensembles_of_point(n_selected_seed)
        integer(int32) :: low_conf_of_point(n_selected_seed)
        integer(int32) :: seed_of_point(1)
        integer(int32) :: n_e, n_se, n_lc, n_seed_of
        integer(int32) :: i, s, r, g, group_id, unit

        call set_ok(ierr)
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_selected_seed, ierr, arg_pos=3_int32, min=0_int32)
        call validate_in_range_int(max_group_size, ierr, arg_pos=4_int32, min=0_int32)
        call validate_in_range_int(n_super_ensembles, ierr, arg_pos=5_int32, min=0_int32)
        if (is_err(ierr)) return

        ensemble_super_id = 0
        do g = 1, n_super_ensembles
            do r = 1, max_group_size
                if (super_ensembles(r, g) == 0) cycle
                ensemble_super_id(super_ensembles(r, g)) = g
            end do
        end do

        vector_to_seed_index = 0
        s = 0
        do i = 1, n_vectors
            if (seed_selection_mask(i)) then
                s = s + 1
                vector_to_seed_index(i) = s
            end if
        end do

        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return

        write (unit, "(A)", advance="no") "row,ensembles,super_ensembles,low_confidence_ensembles,seed_of"//achar(10)

        do i = 1, n_vectors
            n_e = 0
            n_se = 0
            n_lc = 0
            do s = 1, n_selected_seed
                if (ensemble_masks(i, s)) then
                    n_e = n_e + 1
                    ensembles_of_point(n_e) = s

                    group_id = ensemble_super_id(s)
                    if (group_id > 0) then
                        if (stc_int_find(super_ensembles_of_point, n_se, group_id) == 0) then
                            n_se = n_se + 1
                            super_ensembles_of_point(n_se) = group_id
                        end if
                    end if
                end if
                if (ensemble_low_confidence_masks(i, s)) then
                    n_lc = n_lc + 1
                    low_conf_of_point(n_lc) = s
                end if
            end do

            n_seed_of = 0
            if (seed_selection_mask(i)) then
                n_seed_of = 1
                seed_of_point(1) = vector_to_seed_index(i)
            end if

            write (unit, "(I0,',')", advance="no") i
            call stc_write_csv_int_list_quoted(unit, ensembles_of_point, n_e)
            write (unit, "(',')", advance="no")
            call stc_write_csv_int_list_quoted(unit, super_ensembles_of_point, n_se)
            write (unit, "(',')", advance="no")
            call stc_write_csv_int_list_quoted(unit, low_conf_of_point, n_lc)
            write (unit, "(',')", advance="no")
            call stc_write_csv_int_list_quoted(unit, seed_of_point, n_seed_of)
            write (unit, "(A)", advance="no") achar(10)
        end do

        close (unit)
    end subroutine serialize_stc_points_as_csv

    !> M_EXPORT_C
    !| summary: Serializes the full pairwise ensemble Overlap Coefficient matrix as CSV
    !| AUTHOR_ASIS_HALLAB
    !| Only pairs with a nonempty intersection are written, matching `tox_stc_json`'s own
    !| `overlap_coefficient_matrix` convention -- an absent pair means Overlap Coefficient 0.
    subroutine serialize_stc_ensemble_overlap_as_csv(filename, n_vectors, n_selected_seed, ensemble_masks, ierr)
        character(len=*), intent(in) :: filename
            !! Name of the CSV file to write
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        logical(c_bool), intent(in) :: ensemble_masks(n_vectors, n_selected_seed)
            !! Per-ensemble accepted membership, one column per seed
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: i, j, unit, intersect_count

        call set_ok(ierr)
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_selected_seed, ierr, arg_pos=3_int32, min=0_int32)
        if (is_err(ierr)) return

        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return

        write (unit, "(A)", advance="no") "ensemble_a,ensemble_b,overlap_coefficient"//achar(10)

        do i = 1, n_selected_seed - 1
            do j = i + 1, n_selected_seed
                intersect_count = count(ensemble_masks(:, i) .and. ensemble_masks(:, j))
                if (intersect_count < 1) cycle

                write (unit, "(I0,',',I0,',')", advance="no") i, j
                call stc_write_csv_real(unit, real(intersect_count, real64) / &
                                       real(min(count(ensemble_masks(:, i)), count(ensemble_masks(:, j))), real64))
                write (unit, "(A)", advance="no") achar(10)
            end do
        end do

        close (unit)
    end subroutine serialize_stc_ensemble_overlap_as_csv

    !> M_EXPORT_C
    !| summary: Serializes the super-ensembles as a gene-family-file-style TSV
    !| AUTHOR_ASIS_HALLAB
    !| One line per super-ensemble: `<group_id>` TAB `<comma-separated member ensemble ids>`,
    !| no header, no quoting -- see this module's own doc comment for why.
    subroutine serialize_stc_super_ensembles_as_tsv(filename, max_group_size, n_super_ensembles, super_ensembles, ierr)
        character(len=*), intent(in) :: filename
            !! Name of the TSV file to write
        integer(int32), intent(in) :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(int32), intent(in) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        integer(int32), intent(in) :: super_ensembles(max_group_size, n_super_ensembles)
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: g, r, unit, group_size

        call set_ok(ierr)
        call validate_in_range_int(max_group_size, ierr, arg_pos=2_int32, min=0_int32)
        call validate_in_range_int(n_super_ensembles, ierr, arg_pos=3_int32, min=0_int32)
        if (is_err(ierr)) return

        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return

        do g = 1, n_super_ensembles
            group_size = count(super_ensembles(:, g) /= 0)
            write (unit, "(I0,A)", advance="no") g, achar(9)
            do r = 1, group_size
                if (r > 1) write (unit, "(',')", advance="no")
                write (unit, "(I0)", advance="no") super_ensembles(r, g)
            end do
            write (unit, "(A)", advance="no") achar(10)
        end do

        close (unit)
    end subroutine serialize_stc_super_ensembles_as_tsv

    !> First index in `list(1:n)` equal to `val`, or 0 if absent. Small private linear-search
    !| helper -- `n_selected_seed`/`n_super_ensembles` are never large enough (see
    !| `misc/mod_STC.md`'s own complexity notes on this whole family) to justify anything more
    !| than this for a per-point dedup check.
    pure function stc_int_find(list, n, val) result(idx)
        integer(int32), intent(in) :: n, val
        integer(int32), intent(in) :: list(*)
        integer(int32) :: idx
        integer(int32) :: i

        idx = 0
        do i = 1, n
            if (list(i) == val) then
                idx = i
                exit
            end if
        end do
    end function stc_int_find

    !> Writes `list(1:n)` as a double-quoted, comma-separated CSV field (`"1,2,3"`, or `""` if
    !| `n == 0`) -- standard CSV quoting, so a list cell's embedded commas do not collide with
    !| the file's own field separator.
    subroutine stc_write_csv_int_list_quoted(unit, list, n)
        integer(int32), intent(in) :: unit, n
        integer(int32), intent(in) :: list(*)
        integer(int32) :: i

        write (unit, "('""')", advance="no")
        do i = 1, n
            if (i > 1) write (unit, "(',')", advance="no")
            write (unit, "(I0)", advance="no") list(i)
        end do
        write (unit, "('""')", advance="no")
    end subroutine stc_write_csv_int_list_quoted

    !> Writes a real number as a plain CSV field. `NaN`/`Infinite` are written empty (CSV's
    !| own "missing value" convention, parsed as `NaN` by `pandas.read_csv` and `NA` by `read.csv`).
    subroutine stc_write_csv_real(unit, val)
        integer(int32), intent(in) :: unit
        real(real64), intent(in) :: val
        character(len=32) :: buffer

        if (ieee_is_nan(val) .or. .not. ieee_is_finite(val)) return
        write (buffer, "(ES24.16E3)") val
        write (unit, "(A)", advance="no") trim(adjustl(buffer))
    end subroutine stc_write_csv_real

end module tox_stc_csv
