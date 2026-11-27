!> The `get_shift_components` cases: of gene j's shift vector, column j of `shift_vectors`,
!| rows 1..d are the start (its family's centroid) and rows d+1..2d the shift. Refused: a gene
!| index outside 1..n_genes (ERR_INVALID_INPUT), other than 2d rows, and an `out_start` or an
!| `out_shift` of another length than d (ERR_SIZE_MISMATCH), the gene index first.
!|
!| `shift_vectors` must have exactly 2d rows.
module mod_test_get_shift_components
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_negative_inf, ieee_is_nan
    use tox_data_accessors, only: get_shift_components
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: N_SAMPLES, N_GENES, UNWRITTEN, shift_vectors
    implicit none
    private
    public :: get_all_tests_get_shift_components

contains

    !> Every case of `get_shift_components`.
    function get_all_tests_get_shift_components() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_get_shift_components_values", test_get_shift_components_values)
        all_tests(2) = test_case("test_get_shift_components_gene_bound", test_get_shift_components_gene_bound)
        all_tests(3) = test_case("test_get_shift_components_rows", test_get_shift_components_rows)
        all_tests(4) = test_case("test_get_shift_components_lengths", test_get_shift_components_lengths)
        all_tests(5) = test_case("test_get_shift_components_no_sample", test_get_shift_components_no_sample)
        all_tests(6) = test_case("test_get_shift_components_more_rows", test_get_shift_components_more_rows)
        all_tests(7) = test_case("test_get_shift_components_copies_non_finite", &
                                 test_get_shift_components_copies_non_finite)
    end function get_all_tests_get_shift_components

    !> Each gene's start and shift, as the fixture writes them: g1 starts at family 2's
    !| [201, 202] and shifts by [-13, -14], g2 (in no family) at [0, 0] by [-23, -24], g3 at
    !| family 1's [101, 102] by [-33, -34].
    subroutine test_get_shift_components_values()
        real(real64) :: shifts(2*N_SAMPLES, N_GENES), start(N_SAMPLES), shift(N_SAMPLES)
        real(real64) :: expected_starts(N_SAMPLES, N_GENES), expected_shifts(N_SAMPLES, N_GENES)
        integer(int32) :: i_gene, ierr

        shifts = shift_vectors()
        expected_starts = reshape([201.0_real64, 202.0_real64, 0.0_real64, 0.0_real64, 101.0_real64, 102.0_real64], &
                                  [N_SAMPLES, N_GENES])
        expected_shifts = reshape([-13.0_real64, -14.0_real64, -23.0_real64, -24.0_real64, -33.0_real64, -34.0_real64], &
                                  [N_SAMPLES, N_GENES])
        do i_gene = 1, N_GENES
            start = UNWRITTEN
            shift = UNWRITTEN
            call get_shift_components(i_gene, shifts, N_SAMPLES, start, shift, ierr)
            call assert_err(ierr, ERR_OK, "test_get_shift_components_values: ierr of gene "//i_gene)
            call assert_equal_array_real(start, expected_starts(:, i_gene), N_SAMPLES, 0.0_real64, &
                                         "test_get_shift_components_values: start of gene "//i_gene)
            call assert_equal_array_real(shift, expected_shifts(:, i_gene), N_SAMPLES, 0.0_real64, &
                                         "test_get_shift_components_values: shift of gene "//i_gene)
        end do
    end subroutine test_get_shift_components_values

    !> 0 and n_genes + 1 are refused with ERR_INVALID_INPUT and leave both outputs alone; the
    !| gene index is checked before the lengths.
    subroutine test_get_shift_components_gene_bound()
        real(real64) :: shifts(2*N_SAMPLES, N_GENES), start(N_SAMPLES), shift(N_SAMPLES), expected(N_SAMPLES)
        real(real64) :: too_long(N_SAMPLES + 1)
        integer(int32) :: ierr

        shifts = shift_vectors()
        expected = UNWRITTEN
        start = UNWRITTEN
        shift = UNWRITTEN
        call get_shift_components(0, shifts, N_SAMPLES, start, shift, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_shift_components_gene_bound: 0")
        call get_shift_components(N_GENES + 1, shifts, N_SAMPLES, start, shift, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_shift_components_gene_bound: n_genes + 1")
        call assert_equal_array_real(start, expected, N_SAMPLES, 0.0_real64, "test_get_shift_components_gene_bound: start")
        call assert_equal_array_real(shift, expected, N_SAMPLES, 0.0_real64, "test_get_shift_components_gene_bound: shift")
        call get_shift_components(0, shifts, N_SAMPLES, too_long, shift, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_shift_components_gene_bound: 0, and a start too long")
    end subroutine test_get_shift_components_gene_bound

    !> `shift_vectors` needs 2d rows: with d = 2, three rows are refused with ERR_SIZE_MISMATCH,
    !| and so is a d of 3 for the fixture's four rows.
    subroutine test_get_shift_components_rows()
        real(real64) :: shifts(2*N_SAMPLES, N_GENES), short_shifts(2*N_SAMPLES - 1, N_GENES)
        real(real64) :: start(N_SAMPLES), shift(N_SAMPLES), start3(3), shift3(3)
        integer(int32) :: ierr

        shifts = shift_vectors()
        short_shifts = shifts(1:2*N_SAMPLES - 1, :)
        call get_shift_components(1, short_shifts, N_SAMPLES, start, shift, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_rows: 2d - 1 rows")
        call get_shift_components(1, shifts, 3, start3, shift3, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_rows: d = 3 of 4 rows")
    end subroutine test_get_shift_components_rows

    !> `out_start` and `out_shift` each need d values: one more or one fewer in either is refused
    !| with ERR_SIZE_MISMATCH.
    subroutine test_get_shift_components_lengths()
        real(real64) :: shifts(2*N_SAMPLES, N_GENES), right(N_SAMPLES), too_long(N_SAMPLES + 1), too_short(N_SAMPLES - 1)
        integer(int32) :: ierr

        shifts = shift_vectors()
        call get_shift_components(1, shifts, N_SAMPLES, too_long, right, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_lengths: start too long")
        call get_shift_components(1, shifts, N_SAMPLES, too_short, right, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_lengths: start too short")
        call get_shift_components(1, shifts, N_SAMPLES, right, too_long, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_lengths: shift too long")
        call get_shift_components(1, shifts, N_SAMPLES, right, too_short, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_lengths: shift too short")
    end subroutine test_get_shift_components_lengths

    !> No sample, d = 0: a start and a shift of no value, from shift vectors of no row.
    subroutine test_get_shift_components_no_sample()
        real(real64) :: shifts(0, N_GENES), start(0), shift(0)
        integer(int32) :: ierr

        call get_shift_components(1, shifts, 0, start, shift, ierr)
        call assert_err(ierr, ERR_OK, "test_get_shift_components_no_sample")
    end subroutine test_get_shift_components_no_sample

    !> `shift_vectors` has exactly 2d rows: one more is refused with ERR_SIZE_MISMATCH as one fewer
    !| is. (A taller one used to be accepted, its rows past 2d ignored.)
    subroutine test_get_shift_components_more_rows()
        real(real64) :: shifts(2*N_SAMPLES + 1, N_GENES), start(N_SAMPLES), shift(N_SAMPLES)
        integer(int32) :: ierr

        shifts(1:2*N_SAMPLES, :) = shift_vectors()
        shifts(2*N_SAMPLES + 1, :) = 999.0_real64
        call get_shift_components(1, shifts, N_SAMPLES, start, shift, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_shift_components_more_rows")
    end subroutine test_get_shift_components_more_rows

    !> Pins today's behaviour, an open question for the maintainer: the accessors copy what is
    !| stored, NaN and Inf included, with ERR_OK (TOX refuses NaN and Inf as analysis input; the
    !| accessors do not look).
    subroutine test_get_shift_components_copies_non_finite()
        real(real64) :: shifts(2*N_SAMPLES, N_GENES), start(N_SAMPLES), shift(N_SAMPLES)
        integer(int32) :: ierr

        shifts = shift_vectors()
        shifts(1, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        shifts(4, 1) = ieee_value(1.0_real64, ieee_negative_inf)
        call get_shift_components(1, shifts, N_SAMPLES, start, shift, ierr)
        call assert_err(ierr, ERR_OK, "test_get_shift_components_copies_non_finite: ierr")
        call assert_true(ieee_is_nan(start(1)), "test_get_shift_components_copies_non_finite: the NaN")
        call assert_true(shift(2) == ieee_value(1.0_real64, ieee_negative_inf), &
                         "test_get_shift_components_copies_non_finite: the -Inf")
    end subroutine test_get_shift_components_copies_non_finite
end module mod_test_get_shift_components
