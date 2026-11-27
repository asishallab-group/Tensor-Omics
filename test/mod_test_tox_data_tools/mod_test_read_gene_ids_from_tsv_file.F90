!> The `read_gene_ids_from_tsv_file` cases: the ids of one column of a tab-separated file, one per
!| data row and in file order, after `n_header_rows` lines; the header and row-count bounds, an
!| empty file, a file of no data row, and the shipped expression file.
!|
!| Each file is written by the case into a scratch directory, and the expected ids are read off
!| its content.
module mod_test_read_gene_ids_from_tsv_file
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_tools, only: read_gene_ids_from_tsv_file
    use tox_errors
    use test_suite, only: test_case
    use scratch_directory
    use mod_test_tox_data_tools_fixtures, only: TAB, LF
    implicit none
    private
    public :: get_all_tests_read_gene_ids_from_tsv_file

contains

    !> Every case of `read_gene_ids_from_tsv_file`.
    function get_all_tests_read_gene_ids_from_tsv_file() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))
        all_tests(1) = test_case("test_read_gene_ids_from_tsv_file_reads_the_column", &
                                 test_read_gene_ids_from_tsv_file_reads_the_column)
        all_tests(2) = test_case("test_read_gene_ids_from_tsv_file_header_rows", &
                                 test_read_gene_ids_from_tsv_file_header_rows)
        all_tests(3) = test_case("test_read_gene_ids_from_tsv_file_row_count", &
                                 test_read_gene_ids_from_tsv_file_row_count)
        all_tests(4) = test_case("test_read_gene_ids_from_tsv_file_no_data_row", &
                                 test_read_gene_ids_from_tsv_file_no_data_row)
        all_tests(5) = test_case("test_read_gene_ids_from_tsv_file_missing_header", &
                                 test_read_gene_ids_from_tsv_file_missing_header)
        all_tests(6) = test_case("test_read_gene_ids_from_tsv_file_splits_at_tabs_only", &
                                 test_read_gene_ids_from_tsv_file_splits_at_tabs_only)
        all_tests(7) = test_case("test_read_gene_ids_from_tsv_file_keeps_duplicates", &
                                 test_read_gene_ids_from_tsv_file_keeps_duplicates)
        all_tests(8) = test_case("test_read_gene_ids_from_tsv_file_shipped_material", &
                                 test_read_gene_ids_from_tsv_file_shipped_material)
        all_tests(9) = test_case("test_read_gene_ids_from_tsv_file_refuses_an_id_too_long", &
                                 test_read_gene_ids_from_tsv_file_refuses_an_id_too_long)
        all_tests(10) = test_case("test_read_gene_ids_from_tsv_file_reads_a_long_line_whole", &
                                  test_read_gene_ids_from_tsv_file_reads_a_long_line_whole)
    end function get_all_tests_read_gene_ids_from_tsv_file

    !> In place of the old `test_read_gene_ids` (which checked only that the first id of the
    !| shipped file is not empty): three data rows, read from the first column and then from the
    !| third, in file order. Blanks around an id are not part of it, and the last line needs no
    !| line ending.
    subroutine test_read_gene_ids_from_tsv_file_reads_the_column()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=8) :: gene_ids(3), expected(3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_column")
        tsv = scratch_file(scratch, "genes.tsv")
        call write_scratch_file(tsv, "id"//TAB//"a"//TAB//"symbol"//LF// &
                                "g2"//TAB//"1.5"//TAB//"  ABC1"//LF// &
                                " g10 "//TAB//"2.5"//TAB//"XYZ"//LF// &
                                "g1"//TAB//"0"//TAB//"Q")

        call read_gene_ids_from_tsv_file(tsv, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_reads_the_column: column 1, ierr")
        expected = ["g2 ", "g10", "g1 "]
        call assert_equal_array_char(gene_ids, expected, 8, 3, &
                                     "test_read_gene_ids_from_tsv_file_reads_the_column: column 1")

        call read_gene_ids_from_tsv_file(tsv, gene_ids, 1, 3, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_reads_the_column: column 3, ierr")
        expected = ["ABC1", "XYZ ", "Q   "]
        call assert_equal_array_char(gene_ids, expected, 8, 3, &
                                     "test_read_gene_ids_from_tsv_file_reads_the_column: column 3")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_reads_the_column

    !> `n_header_rows` lines are skipped, whatever they hold: none, so the first line is an id,
    !| and two, so the first id is on the third line.
    subroutine test_read_gene_ids_from_tsv_file_header_rows()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=8) :: three_ids(3), one_id(1), expected_three(3), expected_one(1)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_header_rows")
        tsv = scratch_file(scratch, "genes.tsv")
        call write_scratch_file(tsv, "h1"//TAB//"x"//LF//"h2"//LF//"g1"//TAB//"1"//LF)

        call read_gene_ids_from_tsv_file(tsv, three_ids, 0, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_header_rows: none, ierr")
        expected_three = ["h1", "h2", "g1"]
        call assert_equal_array_char(three_ids, expected_three, 8, 3, &
                                     "test_read_gene_ids_from_tsv_file_header_rows: none")

        call read_gene_ids_from_tsv_file(tsv, one_id, 2, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_header_rows: two, ierr")
        expected_one = ["g1"]
        call assert_equal_array_char(one_id, expected_one, 8, 1, "test_read_gene_ids_from_tsv_file_header_rows: two")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_header_rows

    !> Both sides of the row-count bound: three data rows fill three ids, and are refused by two
    !| with ERR_SIZE_MISMATCH. Four ids take three rows (the fourth is not written: the procedure
    !| does not say how many rows it read). And no id at all is too few for one row.
    subroutine test_read_gene_ids_from_tsv_file_row_count()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=8) :: three_ids(3), two_ids(2), four_ids(4), no_ids(0), expected(3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_row_count")
        tsv = scratch_file(scratch, "genes.tsv")
        call write_scratch_file(tsv, "id"//LF//"g1"//LF//"g2"//LF//"g3"//LF)
        expected = ["g1", "g2", "g3"]

        call read_gene_ids_from_tsv_file(tsv, three_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_row_count: as many ids as rows, ierr")
        call assert_equal_array_char(three_ids, expected, 8, 3, &
                                     "test_read_gene_ids_from_tsv_file_row_count: as many ids as rows")

        call read_gene_ids_from_tsv_file(tsv, two_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_gene_ids_from_tsv_file_row_count: one id fewer than rows")

        call read_gene_ids_from_tsv_file(tsv, four_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_row_count: one id more than rows, ierr")
        call assert_equal_array_char(four_ids(1:3), expected, 8, 3, &
                                     "test_read_gene_ids_from_tsv_file_row_count: one id more than rows")

        call read_gene_ids_from_tsv_file(tsv, no_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_gene_ids_from_tsv_file_row_count: no id for three rows")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_row_count

    !> A file of nothing but its header has no gene: read without an error, into no id or into
    !| some. So has an empty file when no header is expected.
    subroutine test_read_gene_ids_from_tsv_file_no_data_row()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: header_only, empty
        character(len=8) :: no_ids(0), two_ids(2)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_no_data_row")
        header_only = scratch_file(scratch, "header_only.tsv")
        call write_scratch_file(header_only, "id"//TAB//"a"//LF)
        empty = scratch_file(scratch, "empty.tsv")
        call write_scratch_file(empty, "")

        call read_gene_ids_from_tsv_file(header_only, no_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_no_data_row: header only, no id")
        call read_gene_ids_from_tsv_file(header_only, two_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_no_data_row: header only, two ids")
        call read_gene_ids_from_tsv_file(empty, no_ids, 0, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_no_data_row: empty file, no header")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_no_data_row

    !> In place of the old empty-file check of `test_read_expression_data`, for this reader: a file
    !| with fewer lines than `n_header_rows` is refused with ERR_READ_DATA. An empty file has no
    !| header line; a one-line file has one header line but not two.
    subroutine test_read_gene_ids_from_tsv_file_missing_header()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: empty, one_line
        character(len=8) :: gene_ids(2)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_missing_header")
        empty = scratch_file(scratch, "empty.tsv")
        call write_scratch_file(empty, "")
        one_line = scratch_file(scratch, "one_line.tsv")
        call write_scratch_file(one_line, "id"//LF)

        call read_gene_ids_from_tsv_file(empty, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_READ_DATA, "test_read_gene_ids_from_tsv_file_missing_header: empty file")
        call read_gene_ids_from_tsv_file(one_line, gene_ids, 2, 1, ierr)
        call assert_err(ierr, ERR_READ_DATA, "test_read_gene_ids_from_tsv_file_missing_header: one line, two expected")
        call read_gene_ids_from_tsv_file(one_line, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_missing_header: one line, one expected")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_missing_header

    !> The reader takes no delimiter: only a tab separates fields, so in a comma-separated file the
    !| whole line is the first field, and the id.
    subroutine test_read_gene_ids_from_tsv_file_splits_at_tabs_only()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: csv
        character(len=8) :: gene_ids(1), expected(1)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_tabs_only")
        csv = scratch_file(scratch, "genes.csv")
        call write_scratch_file(csv, "id,a"//LF//"g1,1.5"//LF)

        call read_gene_ids_from_tsv_file(csv, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_splits_at_tabs_only: ierr")
        expected = ["g1,1.5"]
        call assert_equal_array_char(gene_ids, expected, 8, 1, "test_read_gene_ids_from_tsv_file_splits_at_tabs_only")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_splits_at_tabs_only

    !> Kept from the old `test_read_expression_data` (its file of duplicate ids read without an
    !| error): the reader takes the ids as they are, an id listed twice included. Refusing
    !| duplicates is `validate_string_array_uniqueness`'s job (tox_data_validation).
    subroutine test_read_gene_ids_from_tsv_file_keeps_duplicates()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=8) :: gene_ids(3), expected(3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_duplicates")
        tsv = scratch_file(scratch, "genes.tsv")
        call write_scratch_file(tsv, "id"//LF//"g1"//LF//"g2"//LF//"g1"//LF)

        call read_gene_ids_from_tsv_file(tsv, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_keeps_duplicates: ierr")
        expected = ["g1", "g2", "g1"]
        call assert_equal_array_char(gene_ids, expected, 8, 3, "test_read_gene_ids_from_tsv_file_keeps_duplicates")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_keeps_duplicates

    !> Kept from the old `test_read_gene_ids` and the R suite's `test_1_read_expression_vectors_tsv`
    !| (which checked the count 88327): the shipped expression file has 88327 data rows, so it
    !| fills 88327 ids and is refused by 88326. The first two and the last ids, read off the file
    !| with `awk -F'\t' 'NR==2 || NR==3 || END'`, are NP_001000000.1, NP_001000001.1 and
    !| XP_987278.1.
    subroutine test_read_gene_ids_from_tsv_file_shipped_material()
        character(len=*), parameter :: EXPRESSION_FILE = "material/kallisto_sex_data_no_na.tsv"
        integer(int32), parameter :: N_ROWS = 88327
        character(len=32), allocatable :: gene_ids(:)
        integer(int32) :: ierr

        allocate (gene_ids(N_ROWS))
        call read_gene_ids_from_tsv_file(EXPRESSION_FILE, gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_gene_ids_from_tsv_file_shipped_material: ierr")
        call assert_string_equal(gene_ids(1), "NP_001000000.1", "test_read_gene_ids_from_tsv_file_shipped_material: first")
        call assert_string_equal(gene_ids(2), "NP_001000001.1", "test_read_gene_ids_from_tsv_file_shipped_material: second")
        call assert_string_equal(gene_ids(N_ROWS), "XP_987278.1", "test_read_gene_ids_from_tsv_file_shipped_material: last")

        call read_gene_ids_from_tsv_file(EXPRESSION_FILE, gene_ids(1:N_ROWS - 1), 1, 1, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_gene_ids_from_tsv_file_shipped_material: one id too few")
    end subroutine test_read_gene_ids_from_tsv_file_shipped_material

    !> BUG (fails today): an id longer than the strings of `gene_ids` is cut to fit, without a
    !| word. With versioned ids that is another gene's id: NP_001.10 cut to 8 characters is
    !| NP_001.1. The reader has to either keep the id or refuse it (ERR_STRING_TOO_LONG would say
    !| why); returning the cut id with ERR_OK is wrong.
    subroutine test_read_gene_ids_from_tsv_file_refuses_an_id_too_long()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=8) :: gene_ids(1)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_too_long")
        tsv = scratch_file(scratch, "genes.tsv")
        call write_scratch_file(tsv, "id"//LF//"NP_001.10"//LF)

        call read_gene_ids_from_tsv_file(tsv, gene_ids, 1, 1, ierr)
        call assert_false(is_ok(ierr) .and. gene_ids(1) == "NP_001.1", &
                          "test_read_gene_ids_from_tsv_file_refuses_an_id_too_long: NP_001.10 came back as NP_001.1 with ERR_OK")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_refuses_an_id_too_long

    !> BUG (fails today): a line is read into a buffer of 2048 characters, and what is past it is
    !| dropped without a word. Here the id is in the second column and runs from character 2041 to
    !| 2052: it comes back cut after its eighth character. The reader has to read the whole line
    !| or refuse it; returning the cut id with ERR_OK is wrong.
    subroutine test_read_gene_ids_from_tsv_file_reads_a_long_line_whole()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=16) :: gene_ids(1)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, "gene_ids_long_line")
        tsv = scratch_file(scratch, "genes.tsv")
        ! 2039 characters of padding, the tab at 2040, the id at 2041..2052
        call write_scratch_file(tsv, "h"//TAB//"id"//LF//repeat("x", 2039)//TAB//"GENE00000042"//LF)

        call read_gene_ids_from_tsv_file(tsv, gene_ids, 1, 2, ierr)
        call assert_false(is_ok(ierr) .and. gene_ids(1) /= "GENE00000042", &
                          "test_read_gene_ids_from_tsv_file_reads_a_long_line_whole: the id past character 2048 came back as '" &
                          //trim(gene_ids(1))//"' with ERR_OK")

        call remove_scratch_directory(scratch)
    end subroutine test_read_gene_ids_from_tsv_file_reads_a_long_line_whole
end module mod_test_read_gene_ids_from_tsv_file
