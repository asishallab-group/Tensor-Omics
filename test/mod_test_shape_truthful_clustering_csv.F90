!> Unit test suite for tox_stc_csv (serialize_stc_points_as_csv,
!| serialize_stc_ensemble_overlap_as_csv, serialize_stc_super_ensembles_as_tsv). None of these
!| three writers touch the observable/history arrays -- the fixture only needs
!| seed_selection_mask/ensemble_masks/ensemble_low_confidence_masks/super_ensembles, the same
!| ones test_json_two_ensembles_with_overlap in
!| test/mod_test_shape_truthful_clustering_json.F90 uses for its own membership assertions.
module mod_test_shape_truthful_clustering_csv
    use tox_stc_csv, only: serialize_stc_points_as_csv, serialize_stc_ensemble_overlap_as_csv, &
        serialize_stc_super_ensembles_as_tsv
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_csv() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(6))

        all_tests(1) = test_case("test_points_csv_membership", test_points_csv_membership)
        all_tests(2) = test_case("test_points_csv_zero_ensembles", test_points_csv_zero_ensembles)
        all_tests(3) = test_case("test_overlap_csv_pairwise_coefficient", test_overlap_csv_pairwise_coefficient)
        all_tests(4) = test_case("test_overlap_csv_no_intersecting_pairs", test_overlap_csv_no_intersecting_pairs)
        all_tests(5) = test_case("test_super_ensembles_tsv_gene_family_format", &
                                 test_super_ensembles_tsv_gene_family_format)
        all_tests(6) = test_case("test_super_ensembles_tsv_empty_when_zero_groups", &
                                 test_super_ensembles_tsv_empty_when_zero_groups)
    end function get_all_tests_shape_truthful_clustering_csv

    !> N=4, 2 ensembles: {1,2,3} (seed=1) and {2,3,4} (seed=4), overlapping on {2,3} -- Overlap
    !| Coefficient 2/3 -- merged into one super-ensemble. Point 1 is also flagged
    !| low-confidence. Matches
    !| mod_test_shape_truthful_clustering_json.F90:test_json_two_ensembles_with_overlap's own
    !| membership fixture exactly.
    subroutine build_fixture(seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles)
        logical(c_bool), intent(out) :: seed_selection_mask(4)
        logical(c_bool), intent(out) :: ensemble_masks(4, 2)
        logical(c_bool), intent(out) :: ensemble_low_confidence_masks(4, 2)
        integer(int32), intent(out) :: super_ensembles(2, 2)

        seed_selection_mask = [.true., .false., .false., .true.]
        ensemble_masks(:, 1) = [.true., .true., .true., .false.]
        ensemble_masks(:, 2) = [.false., .true., .true., .true.]
        ensemble_low_confidence_masks = .false.
        ensemble_low_confidence_masks(1, 1) = .true.
        super_ensembles(:, 1) = [1, 2]
        super_ensembles(:, 2) = [0, 0]
    end subroutine build_fixture

    subroutine test_points_csv_membership()
        logical(c_bool) :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        integer(int32) :: super_ensembles(2, 2)
        character(len=32) :: filename
        character(len=200) :: line
        integer(int32) :: ierr, unit

        call build_fixture(seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles)

        filename = 'test_stc_points.csv'
        call serialize_stc_points_as_csv(filename, 4_int32, 2_int32, 2_int32, 1_int32, &
                                         seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, &
                                         super_ensembles, ierr)
        call assert_true(is_ok(ierr), "points csv: must not fail")

        open (newunit=unit, file=trim(filename), status='old', action='read')

        read (unit, "(A)") line
        call assert_string_equal(trim(line), "row,ensembles,super_ensembles,low_confidence_ensembles,seed_of", &
                                 "points csv: header")

        read (unit, "(A)") line
        call assert_string_equal(trim(line), '1,"1","1","1","1"', "points csv: row 1 (seed, low-confidence)")

        read (unit, "(A)") line
        call assert_string_equal(trim(line), '2,"1,2","1","",""', "points csv: row 2 (in both ensembles)")

        read (unit, "(A)") line
        call assert_string_equal(trim(line), '3,"1,2","1","",""', "points csv: row 3 (in both ensembles)")

        read (unit, "(A)") line
        call assert_string_equal(trim(line), '4,"2","1","","2"', "points csv: row 4 (seed of ensemble 2)")

        close (unit, status='delete')
    end subroutine test_points_csv_membership

    subroutine test_points_csv_zero_ensembles()
        logical(c_bool) :: seed_selection_mask(2), ensemble_masks(2, 0), ensemble_low_confidence_masks(2, 0)
        integer(int32) :: super_ensembles(2, 0)
        character(len=32) :: filename
        character(len=200) :: line
        integer(int32) :: ierr, unit

        seed_selection_mask = .false.

        filename = 'test_stc_points_zero.csv'
        call serialize_stc_points_as_csv(filename, 2_int32, 0_int32, 2_int32, 0_int32, &
                                         seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, &
                                         super_ensembles, ierr)
        call assert_true(is_ok(ierr), "points csv zero ensembles: must not fail")

        open (newunit=unit, file=trim(filename), status='old', action='read')
        read (unit, "(A)") line
        read (unit, "(A)") line
        call assert_string_equal(trim(line), '1,"","","",""', "points csv zero ensembles: all-empty row")
        close (unit, status='delete')
    end subroutine test_points_csv_zero_ensembles

    subroutine test_overlap_csv_pairwise_coefficient()
        logical(c_bool) :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        integer(int32) :: super_ensembles(2, 2)
        character(len=32) :: filename
        character(len=200) :: line
        integer(int32) :: ierr, unit

        call build_fixture(seed_selection_mask, ensemble_masks, ensemble_low_confidence_masks, super_ensembles)

        filename = 'test_stc_overlap.csv'
        call serialize_stc_ensemble_overlap_as_csv(filename, 4_int32, 2_int32, ensemble_masks, ierr)
        call assert_true(is_ok(ierr), "overlap csv: must not fail")

        open (newunit=unit, file=trim(filename), status='old', action='read')
        read (unit, "(A)") line
        call assert_string_equal(trim(line), "ensemble_a,ensemble_b,overlap_coefficient", "overlap csv: header")
        read (unit, "(A)") line
        call assert_string_equal(trim(line), "1,2,6.6666666666666663E-001", "overlap csv: the one intersecting pair")
        close (unit, status='delete')
    end subroutine test_overlap_csv_pairwise_coefficient

    subroutine test_overlap_csv_no_intersecting_pairs()
        logical(c_bool) :: ensemble_masks(4, 2)
        character(len=32) :: filename
        character(len=200) :: line
        integer(int32) :: ierr, unit, ios

        ensemble_masks(:, 1) = [.true., .true., .false., .false.]
        ensemble_masks(:, 2) = [.false., .false., .true., .true.]

        filename = 'test_stc_overlap_disjoint.csv'
        call serialize_stc_ensemble_overlap_as_csv(filename, 4_int32, 2_int32, ensemble_masks, ierr)
        call assert_true(is_ok(ierr), "overlap csv disjoint: must not fail")

        open (newunit=unit, file=trim(filename), status='old', action='read')
        read (unit, "(A)") line
        call assert_string_equal(trim(line), "ensemble_a,ensemble_b,overlap_coefficient", "overlap csv disjoint: header")
        ! either straight EOF, or one blank record from gfortran's own trailing-newline-on-close
        ! quirk (see tox_stc_html_assets's doc comment) -- either way, no actual data row
        line = ''
        read (unit, "(A)", iostat=ios) line
        call assert_true(ios /= 0 .or. len_trim(line) == 0, &
                         "overlap csv disjoint: no data rows, since no pair intersects")
        close (unit, status='delete')
    end subroutine test_overlap_csv_no_intersecting_pairs

    subroutine test_super_ensembles_tsv_gene_family_format()
        integer(int32) :: super_ensembles(3, 1)
        character(len=32) :: filename
        character(len=200) :: line
        integer(int32) :: ierr, unit

        super_ensembles(:, 1) = [1, 2, 0]

        filename = 'test_stc_super_ensembles.tsv'
        call serialize_stc_super_ensembles_as_tsv(filename, 3_int32, 1_int32, super_ensembles, ierr)
        call assert_true(is_ok(ierr), "super_ensembles tsv: must not fail")

        open (newunit=unit, file=trim(filename), status='old', action='read')
        read (unit, "(A)") line
        call assert_string_equal(trim(line), "1"//achar(9)//"1,2", "super_ensembles tsv: id TAB comma-list, no quoting")
        close (unit, status='delete')
    end subroutine test_super_ensembles_tsv_gene_family_format

    subroutine test_super_ensembles_tsv_empty_when_zero_groups()
        integer(int32) :: super_ensembles(2, 0)
        character(len=32) :: filename
        integer(int32) :: ierr, filesize, unit

        filename = 'test_stc_super_ensembles_empty.tsv'
        call serialize_stc_super_ensembles_as_tsv(filename, 2_int32, 0_int32, super_ensembles, ierr)
        call assert_true(is_ok(ierr), "super_ensembles tsv empty: must not fail")

        inquire (file=trim(filename), size=filesize)
        ! at most the harmless single trailing-newline-on-close byte, no group lines
        call assert_true(filesize <= 1, "super_ensembles tsv empty: no group lines written")
        open (newunit=unit, file=trim(filename), status='old')
        close (unit, status='delete')
    end subroutine test_super_ensembles_tsv_empty_when_zero_groups

end module mod_test_shape_truthful_clustering_csv
