!> The `read_expression_vectors_tsv` cases: each data row's values, taken from `value_cols`, go to
!| the column of `expression_vectors` whose gene the row names, starting at row `start_row`, and
!| each further file continues below the last; the delimiter, the header and every refusal the
!| reader makes (no value column, the gene count, a missing header, a file of no data row, a
!| value column past the file's columns, a short row, a value that is not a number or not finite,
!| more rows than genes).
!|
!| Each file is written by the case into a scratch directory, and the expected values are read
!| off its content. They are decimals that binary represents exactly, so every comparison has
!| tolerance 0. `expression_vectors` starts at UNWRITTEN, so a value the reader leaves alone
!| shows.
module mod_test_read_expression_vectors_tsv
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_data_tools, only: read_expression_vectors_tsv
    use tox_errors
    use test_suite, only: test_case
    use scratch_directory
    use mod_test_tox_data_tools_fixtures, only: TAB, LF, UNWRITTEN
    implicit none
    private
    public :: get_all_tests_read_expression_vectors_tsv

    !> The genes most cases read.
    character(len=4), parameter :: THREE_GENE_IDS(3) = ["g1", "g2", "g3"]

contains

    !> Every case of `read_expression_vectors_tsv`.
    function get_all_tests_read_expression_vectors_tsv() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(27))
        all_tests(1) = test_case("test_read_expression_vectors_tsv_values", test_read_expression_vectors_tsv_values)
        all_tests(2) = test_case("test_read_expression_vectors_tsv_value_cols_order", &
                                 test_read_expression_vectors_tsv_value_cols_order)
        all_tests(3) = test_case("test_read_expression_vectors_tsv_start_row", test_read_expression_vectors_tsv_start_row)
        all_tests(4) = test_case("test_read_expression_vectors_tsv_several_files", &
                                 test_read_expression_vectors_tsv_several_files)
        all_tests(5) = test_case("test_read_expression_vectors_tsv_gene_col", test_read_expression_vectors_tsv_gene_col)
        all_tests(6) = test_case("test_read_expression_vectors_tsv_unmatched_genes", &
                                 test_read_expression_vectors_tsv_unmatched_genes)
        all_tests(7) = test_case("test_read_expression_vectors_tsv_delimiter", test_read_expression_vectors_tsv_delimiter)
        all_tests(8) = test_case("test_read_expression_vectors_tsv_header_rows", test_read_expression_vectors_tsv_header_rows)
        all_tests(9) = test_case("test_read_expression_vectors_tsv_extreme_values", &
                                 test_read_expression_vectors_tsv_extreme_values)
        all_tests(10) = test_case("test_read_expression_vectors_tsv_no_file", test_read_expression_vectors_tsv_no_file)
        all_tests(11) = test_case("test_read_expression_vectors_tsv_refuses_no_value_column", &
                                  test_read_expression_vectors_tsv_refuses_no_value_column)
        all_tests(12) = test_case("test_read_expression_vectors_tsv_gene_count", test_read_expression_vectors_tsv_gene_count)
        all_tests(13) = test_case("test_read_expression_vectors_tsv_missing_header", &
                                  test_read_expression_vectors_tsv_missing_header)
        all_tests(14) = test_case("test_read_expression_vectors_tsv_no_data_row", &
                                  test_read_expression_vectors_tsv_no_data_row)
        all_tests(15) = test_case("test_read_expression_vectors_tsv_value_column_bound", &
                                  test_read_expression_vectors_tsv_value_column_bound)
        all_tests(16) = test_case("test_read_expression_vectors_tsv_refuses_a_short_row", &
                                  test_read_expression_vectors_tsv_refuses_a_short_row)
        all_tests(17) = test_case("test_read_expression_vectors_tsv_refuses_non_numbers", &
                                  test_read_expression_vectors_tsv_refuses_non_numbers)
        all_tests(18) = test_case("test_read_expression_vectors_tsv_refuses_nan_and_inf", &
                                  test_read_expression_vectors_tsv_refuses_nan_and_inf)
        all_tests(19) = test_case("test_read_expression_vectors_tsv_row_count", test_read_expression_vectors_tsv_row_count)
        all_tests(20) = test_case("test_read_expression_vectors_tsv_trailing_empty_line", &
                                  test_read_expression_vectors_tsv_trailing_empty_line)
        all_tests(21) = test_case("test_read_expression_vectors_tsv_refuses_duplicate_rows", &
                                  test_read_expression_vectors_tsv_refuses_duplicate_rows)
        all_tests(22) = test_case("test_read_expression_vectors_tsv_refuses_an_empty_cell", &
                                  test_read_expression_vectors_tsv_refuses_an_empty_cell)
        all_tests(23) = test_case("test_read_expression_vectors_tsv_refuses_a_partial_number", &
                                  test_read_expression_vectors_tsv_refuses_a_partial_number)
        all_tests(24) = test_case("test_read_expression_vectors_tsv_matches_whole_ids", &
                                  test_read_expression_vectors_tsv_matches_whole_ids)
        all_tests(25) = test_case("test_read_expression_vectors_tsv_reads_a_long_line_whole", &
                                  test_read_expression_vectors_tsv_reads_a_long_line_whole)
        all_tests(26) = test_case("test_read_expression_vectors_tsv_missing_file", &
                                  test_read_expression_vectors_tsv_missing_file)
        all_tests(27) = test_case("test_read_expression_vectors_tsv_checks_the_positions", &
                                  test_read_expression_vectors_tsv_checks_the_positions)
    end function get_all_tests_read_expression_vectors_tsv

    !> Writes `content` to `name` in `scratch`, reads it as the one file with one header row and
    !| the genes in column 1, into `expression_vectors`, which starts at UNWRITTEN, and checks
    !| `ierr` against `expected_code`.
    subroutine read_one_file(scratch, name, content, gene_ids, value_cols, expression_vectors, expected_code, msg)
        type(scratch_directory_type), intent(inout) :: scratch
        character(*), intent(in) :: name, content, gene_ids(:), msg
        integer(int32), intent(in) :: value_cols(:), expected_code
        real(real64), intent(inout) :: expression_vectors(:, :)
        character(len=:), allocatable :: tsv
        character(len=1024) :: file_list(1)
        integer(int32) :: ierr

        tsv = scratch_file(scratch, name)
        call write_scratch_file(tsv, content)
        file_list(1) = tsv
        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, expected_code, msg//": ierr")
    end subroutine read_one_file

    !> In place of the old `test_read_expression_data`'s checks of the shipped file (its shape,
    !| and that no value is negative): three rows naming the genes in another order than
    !| `gene_ids`, so each row's values land in its gene's column, not in its row's. Columns 2
    !| and 3 are read, column 4 is not. The reader does not judge values: a negative one is read
    !| as it is.
    subroutine test_read_expression_vectors_tsv_values()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(2, 3), expected(2, 3)
        integer(int32) :: value_cols(2)

        call create_scratch_directory(scratch, "expression_values")
        value_cols = [2, 3]
        call read_one_file(scratch, "expression.tsv", &
                           "id"//TAB//"a"//TAB//"b"//TAB//"unused"//LF// &
                           "g2"//TAB//"1.5"//TAB//"-2"//TAB//"9"//LF// &
                           "g3"//TAB//"0.25"//TAB//"4"//TAB//"9"//LF// &
                           "g1"//TAB//"8"//TAB//"0.5"//TAB//"9"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_OK, "test_read_expression_vectors_tsv_values")
        expected = reshape([8.0_real64, 0.5_real64, 1.5_real64, -2.0_real64, 0.25_real64, 4.0_real64], [2, 3])
        call assert_equal_array_real(expression_vectors, expected, 6, 0.0_real64, "test_read_expression_vectors_tsv_values", 2)

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_values

    !> Row k of a gene's vector is the k-th entry of `value_cols`, in that order: [4, 2] puts
    !| column 4 in row 1 and column 2 in row 2, and [3] reads column 3 alone.
    subroutine test_read_expression_vectors_tsv_value_cols_order()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: content
        real(real64) :: two_rows(2, 3), one_row(1, 3), expected_two(2, 3), expected_one(1, 3)
        integer(int32) :: reversed(2), third(1)

        call create_scratch_directory(scratch, "expression_value_cols_order")
        content = "id"//TAB//"a"//TAB//"b"//TAB//"c"//LF// &
                  "g1"//TAB//"1"//TAB//"2"//TAB//"3"//LF// &
                  "g2"//TAB//"4"//TAB//"5"//TAB//"6"//LF// &
                  "g3"//TAB//"7"//TAB//"8"//TAB//"9"//LF
        reversed = [4, 2]
        call read_one_file(scratch, "expression.tsv", content, THREE_GENE_IDS, reversed, two_rows, ERR_OK, &
                           "test_read_expression_vectors_tsv_value_cols_order: [4, 2]")
        expected_two = reshape([3.0_real64, 1.0_real64, 6.0_real64, 4.0_real64, 9.0_real64, 7.0_real64], [2, 3])
        call assert_equal_array_real(two_rows, expected_two, 6, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_value_cols_order: [4, 2]", 2)

        third = [3]
        call read_one_file(scratch, "expression.tsv", content, THREE_GENE_IDS, third, one_row, ERR_OK, &
                           "test_read_expression_vectors_tsv_value_cols_order: [3]")
        expected_one = reshape([2.0_real64, 5.0_real64, 8.0_real64], [1, 3])
        call assert_equal_array_real(one_row, expected_one, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_value_cols_order: [3]")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_value_cols_order

    !> `start_row` 2 writes the two value columns to rows 2 and 3, and leaves row 1 alone.
    subroutine test_read_expression_vectors_tsv_start_row()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=1024) :: file_list(1)
        real(real64) :: expression_vectors(3, 3), expected(3, 3)
        integer(int32) :: value_cols(2), ierr

        call create_scratch_directory(scratch, "expression_start_row")
        tsv = scratch_file(scratch, "expression.tsv")
        call write_scratch_file(tsv, "id"//TAB//"a"//TAB//"b"//LF// &
                                "g1"//TAB//"1"//TAB//"2"//LF//"g2"//TAB//"3"//TAB//"4"//LF//"g3"//TAB//"5"//TAB//"6"//LF)
        file_list(1) = tsv
        value_cols = [2, 3]
        expression_vectors = UNWRITTEN

        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 2, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_start_row: ierr")
        expected = reshape([UNWRITTEN, 1.0_real64, 2.0_real64, UNWRITTEN, 3.0_real64, 4.0_real64, &
                            UNWRITTEN, 5.0_real64, 6.0_real64], [3, 3])
        call assert_equal_array_real(expression_vectors, expected, 9, 0.0_real64, "test_read_expression_vectors_tsv_start_row", 3)

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_start_row

    !> Two files, as the shipped data comes per tissue: the first file's two value columns fill
    !| rows 1 and 2, the second file's rows 3 and 4. Each file names the genes in its own order,
    !| and the first lacks g3, whose rows 1 and 2 stay as they were.
    subroutine test_read_expression_vectors_tsv_several_files()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: file_list(2)
        real(real64) :: expression_vectors(4, 3), expected(4, 3)
        integer(int32) :: value_cols(2), ierr

        call create_scratch_directory(scratch, "expression_several_files")
        file_list(1) = scratch_file(scratch, "adipose.tsv")
        call write_scratch_file(trim(file_list(1)), "id"//TAB//"r1"//TAB//"r2"//LF// &
                                "g2"//TAB//"1"//TAB//"2"//LF//"g1"//TAB//"3"//TAB//"4"//LF)
        file_list(2) = scratch_file(scratch, "brain.tsv")
        call write_scratch_file(trim(file_list(2)), "id"//TAB//"r1"//TAB//"r2"//LF// &
                                "g3"//TAB//"5"//TAB//"8"//LF//"g1"//TAB//"6"//TAB//"9"//LF//"g2"//TAB//"7"//TAB//"10"//LF)
        value_cols = [2, 3]
        expression_vectors = UNWRITTEN

        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_several_files: ierr")
        expected = reshape([3.0_real64, 4.0_real64, 6.0_real64, 9.0_real64, &
                            1.0_real64, 2.0_real64, 7.0_real64, 10.0_real64, &
                            UNWRITTEN, UNWRITTEN, 5.0_real64, 8.0_real64], [4, 3])
        call assert_equal_array_real(expression_vectors, expected, 12, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_several_files", 4)

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_several_files

    !> The ids may be in any column: here the third, with the values before them.
    subroutine test_read_expression_vectors_tsv_gene_col()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: file_list(1)
        real(real64) :: expression_vectors(2, 3), expected(2, 3)
        integer(int32) :: value_cols(2), ierr

        call create_scratch_directory(scratch, "expression_gene_col")
        file_list(1) = scratch_file(scratch, "expression.tsv")
        call write_scratch_file(trim(file_list(1)), "a"//TAB//"b"//TAB//"id"//LF// &
                                "1"//TAB//"2"//TAB//"g3"//LF//"3"//TAB//"4"//TAB//"g1"//LF//"5"//TAB//"6"//TAB//"g2"//LF)
        value_cols = [1, 2]
        expression_vectors = UNWRITTEN

        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 3, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_gene_col: ierr")
        expected = reshape([3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64, 1.0_real64, 2.0_real64], [2, 3])
        call assert_equal_array_real(expression_vectors, expected, 6, 0.0_real64, "test_read_expression_vectors_tsv_gene_col", 2)

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_gene_col

    !> Kept from the old `test_read_expression_data` (a file whose rows name no known gene reads
    !| without an error): a row of a gene not in `gene_ids` is skipped, and a gene no row names
    !| keeps what its column held.
    subroutine test_read_expression_vectors_tsv_unmatched_genes()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1)

        call create_scratch_directory(scratch, "expression_unmatched_genes")
        value_cols = [2]
        call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF//"gX"//TAB//"2"//LF// &
                           "g3"//TAB//"3"//LF, THREE_GENE_IDS, value_cols, expression_vectors, ERR_OK, &
                           "test_read_expression_vectors_tsv_unmatched_genes")
        expected = reshape([1.0_real64, UNWRITTEN, 3.0_real64], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_unmatched_genes")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_unmatched_genes

    !> Kept from the old `test_read_expression_data` (the CSV file): with `delimiter` ',' a
    !| comma-separated file is read, blanks after a comma included. Without it the delimiter is a
    !| tab, so the same file has one field per row and is refused with ERR_INVALID_INPUT.
    subroutine test_read_expression_vectors_tsv_delimiter()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: file_list(1)
        real(real64) :: expression_vectors(2, 3), expected(2, 3)
        integer(int32) :: value_cols(2), ierr

        call create_scratch_directory(scratch, "expression_delimiter")
        file_list(1) = scratch_file(scratch, "expression.csv")
        call write_scratch_file(trim(file_list(1)), "id,a,b"//LF//"g1,1.5, 2"//LF//"g2, 3,4"//LF//"g3,5,6.25"//LF)
        value_cols = [2, 3]

        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ",", ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_delimiter: comma, ierr")
        expected = reshape([1.5_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.25_real64], [2, 3])
        call assert_equal_array_real(expression_vectors, expected, 6, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_delimiter: comma", 2)

        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_delimiter: the default, a tab")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_delimiter

    !> `n_header_rows` lines are skipped: none, and two.
    subroutine test_read_expression_vectors_tsv_header_rows()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: no_header(1), two_headers(1)
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1), ierr

        call create_scratch_directory(scratch, "expression_header_rows")
        no_header(1) = scratch_file(scratch, "no_header.tsv")
        call write_scratch_file(trim(no_header(1)), "g1"//TAB//"1"//LF//"g2"//TAB//"2"//LF//"g3"//TAB//"3"//LF)
        two_headers(1) = scratch_file(scratch, "two_headers.tsv")
        call write_scratch_file(trim(two_headers(1)), "id"//TAB//"a"//LF//"unit"//TAB//"TPM"//LF// &
                                "g1"//TAB//"4"//LF//"g2"//TAB//"5"//LF//"g3"//TAB//"6"//LF)
        value_cols = [2]

        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(no_header, THREE_GENE_IDS, expression_vectors, 0, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_header_rows: none, ierr")
        expected = reshape([1.0_real64, 2.0_real64, 3.0_real64], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_header_rows: none")

        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(two_headers, THREE_GENE_IDS, expression_vectors, 2, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_header_rows: two, ierr")
        expected = reshape([4.0_real64, 5.0_real64, 6.0_real64], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_header_rows: two")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_header_rows

    !> The ends of the finite range read exactly: the largest finite value written out to 17
    !| digits is huge(1.0_real64), its negative is -huge, and 2.2250738585072014e-308 is
    !| tiny(1.0_real64), the smallest normal value.
    subroutine test_read_expression_vectors_tsv_extreme_values()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1)

        call create_scratch_directory(scratch, "expression_extreme_values")
        value_cols = [2]
        call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF// &
                           "g1"//TAB//"1.7976931348623157e308"//LF// &
                           "g2"//TAB//"-1.7976931348623157e308"//LF// &
                           "g3"//TAB//"2.2250738585072014e-308"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_OK, "test_read_expression_vectors_tsv_extreme_values")
        expected = reshape([huge(1.0_real64), -huge(1.0_real64), tiny(1.0_real64)], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_extreme_values")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_extreme_values

    !> No file at all reads nothing, without an error.
    subroutine test_read_expression_vectors_tsv_no_file()
        character(len=8) :: file_list(0)
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1), ierr

        value_cols = [2]
        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_no_file: ierr")
        expected = UNWRITTEN
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, "test_read_expression_vectors_tsv_no_file")
    end subroutine test_read_expression_vectors_tsv_no_file

    !> No value column is refused with ERR_INVALID_INPUT, before any file is opened (this one
    !| does not exist).
    subroutine test_read_expression_vectors_tsv_refuses_no_value_column()
        character(len=32) :: file_list(1)
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(0), ierr

        file_list(1) = "no_such_expression_file.tsv"
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_refuses_no_value_column")
    end subroutine test_read_expression_vectors_tsv_refuses_no_value_column

    !> `expression_vectors` has a column per gene: one column more or fewer than three genes is
    !| refused with ERR_SIZE_MISMATCH, before any file is opened (this one does not exist).
    subroutine test_read_expression_vectors_tsv_gene_count()
        character(len=32) :: file_list(1)
        real(real64) :: one_more(1, 4), one_fewer(1, 2)
        integer(int32) :: value_cols(1), ierr

        file_list(1) = "no_such_expression_file.tsv"
        value_cols = [2]
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, one_more, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_expression_vectors_tsv_gene_count: a column more than genes")
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, one_fewer, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_expression_vectors_tsv_gene_count: a column fewer than genes")
    end subroutine test_read_expression_vectors_tsv_gene_count

    !> Kept from the old `test_read_expression_data` (the empty file): a file with fewer lines
    !| than `n_header_rows` is refused with ERR_READ_DATA: an empty one, and one line where two
    !| are expected.
    subroutine test_read_expression_vectors_tsv_missing_header()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: empty(1), one_line(1)
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(1), ierr

        call create_scratch_directory(scratch, "expression_missing_header")
        empty(1) = scratch_file(scratch, "empty.tsv")
        call write_scratch_file(trim(empty(1)), "")
        one_line(1) = scratch_file(scratch, "one_line.tsv")
        call write_scratch_file(trim(one_line(1)), "id"//TAB//"a"//LF)
        value_cols = [2]

        call read_expression_vectors_tsv(empty, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_READ_DATA, "test_read_expression_vectors_tsv_missing_header: empty file")
        call read_expression_vectors_tsv(one_line, THREE_GENE_IDS, expression_vectors, 2, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_READ_DATA, "test_read_expression_vectors_tsv_missing_header: one line, two expected")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_missing_header

    !> A file of no data row is refused with ERR_EMPTY_INPUT naming file_list (argument 1): one of
    !| nothing but its header, one of its header and blank lines, and an empty one read without a
    !| header. (The first used to be ERR_READ_DATA, while the gene-id reader read it as no gene.)
    subroutine test_read_expression_vectors_tsv_no_data_row()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: file_list(1)
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(1), ierr

        call create_scratch_directory(scratch, "expression_no_data_row")
        value_cols = [2]
        call read_one_file(scratch, "header_only.tsv", "id"//TAB//"a"//LF, THREE_GENE_IDS, value_cols, expression_vectors, &
                           ERR_EMPTY_INPUT, "test_read_expression_vectors_tsv_no_data_row: header only")
        call read_one_file(scratch, "header_and_blanks.tsv", "id"//TAB//"a"//LF//LF//"  "//LF, THREE_GENE_IDS, value_cols, &
                           expression_vectors, ERR_EMPTY_INPUT, "test_read_expression_vectors_tsv_no_data_row: header, blanks")
        file_list(1) = scratch_file(scratch, "empty.tsv")
        call write_scratch_file(trim(file_list(1)), "")
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 0, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_read_expression_vectors_tsv_no_data_row: empty, no header", 1)

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_no_data_row

    !> A value column must be one of the file's columns, counted on its first data row: of three
    !| columns, 3 is read and 4 is refused with ERR_INVALID_INPUT. With two files the bound is
    !| each file's own, so column 3 is refused when the second file has two columns.
    subroutine test_read_expression_vectors_tsv_value_column_bound()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: two_files(2)
        real(real64) :: expression_vectors(1, 3), expected(1, 3), two_rows(2, 3)
        integer(int32) :: last(1), past(1), ierr

        call create_scratch_directory(scratch, "expression_value_column_bound")
        last = [3]
        call read_one_file(scratch, "three_columns.tsv", "id"//TAB//"a"//TAB//"b"//LF// &
                           "g1"//TAB//"1"//TAB//"2"//LF//"g2"//TAB//"3"//TAB//"4"//LF//"g3"//TAB//"5"//TAB//"6"//LF, &
                           THREE_GENE_IDS, last, expression_vectors, ERR_OK, "test_read_expression_vectors_tsv_value_column_bound: 3")
        expected = reshape([2.0_real64, 4.0_real64, 6.0_real64], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_value_column_bound: 3")

        past = [4]
        call read_one_file(scratch, "three_columns.tsv", "id"//TAB//"a"//TAB//"b"//LF// &
                           "g1"//TAB//"1"//TAB//"2"//LF//"g2"//TAB//"3"//TAB//"4"//LF//"g3"//TAB//"5"//TAB//"6"//LF, &
                           THREE_GENE_IDS, past, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_value_column_bound: 4")

        two_files(1) = scratch_file(scratch, "three_columns.tsv")
        two_files(2) = scratch_file(scratch, "two_columns.tsv")
        call write_scratch_file(trim(two_files(2)), "id"//TAB//"a"//LF//"g1"//TAB//"7"//LF)
        call read_expression_vectors_tsv(two_files, THREE_GENE_IDS, two_rows, 1, 1, last, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_value_column_bound: 3, second file")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_value_column_bound

    !> Kept from the old `test_read_expression_data` (the mixed-separator file; also its unused
    !| missing-value file): a row with fewer fields than the value columns need is refused with
    !| ERR_INVALID_INPUT, wherever it is. Here the last row separates its last two values with a
    !| comma, so it has one field too few; and a row missing its last value.
    subroutine test_read_expression_vectors_tsv_refuses_a_short_row()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(2, 3)
        integer(int32) :: value_cols(2)

        call create_scratch_directory(scratch, "expression_short_row")
        value_cols = [2, 3]
        call read_one_file(scratch, "mixed_separators.tsv", "id"//TAB//"a"//TAB//"b"//LF// &
                           "g1"//TAB//"1"//TAB//"2"//LF//"g2"//TAB//"3"//TAB//"4"//LF//"g3"//TAB//"5, 6"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_refuses_a_short_row: mixed separators")
        call read_one_file(scratch, "missing_value.tsv", "id"//TAB//"a"//TAB//"b"//LF// &
                           "g1"//TAB//"1"//TAB//"2"//LF//"g2"//TAB//"3"//LF//"g3"//TAB//"5"//TAB//"6"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_refuses_a_short_row: a missing value")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_a_short_row

    !> Kept from the old `test_read_expression_data` (the file whose second column held the ids):
    !| a value that is no number is refused with ERR_READ_DATA. An id, "NA", and a word.
    subroutine test_read_expression_vectors_tsv_refuses_non_numbers()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(1), i_word
        character(len=8), parameter :: WORDS(3) = ["g1   ", "NA   ", "high "]

        call create_scratch_directory(scratch, "expression_non_numbers")
        value_cols = [2]
        do i_word = 1, size(WORDS)
            call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF// &
                               "g2"//TAB//trim(WORDS(i_word))//LF//"g3"//TAB//"3"//LF, &
                               THREE_GENE_IDS, value_cols, expression_vectors, ERR_READ_DATA, &
                               "test_read_expression_vectors_tsv_refuses_non_numbers: "//trim(WORDS(i_word)))
        end do

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_non_numbers

    !> Kept from the old `test_read_expression_data` (the NaN and the Inf files), with every
    !| spelling Fortran reads as a non-finite value: refused with ERR_INVALID_INPUT, as TOX takes
    !| no NaN or Inf. A number past huge is refused too; gfortran reads it as Inf (so
    !| ERR_INVALID_INPUT), ifx fails to read it (ERR_READ_DATA), so that one checks for an error.
    subroutine test_read_expression_vectors_tsv_refuses_nan_and_inf()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: tsv
        character(len=1024) :: file_list(1)
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(1), i_word, ierr
        character(len=9), parameter :: WORDS(6) = ["NaN      ", "nan      ", "Inf      ", "-Inf     ", &
                                                   "Infinity ", "+Infinity"]

        call create_scratch_directory(scratch, "expression_nan_and_inf")
        value_cols = [2]
        do i_word = 1, size(WORDS)
            call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF// &
                               "g2"//TAB//trim(WORDS(i_word))//LF//"g3"//TAB//"3"//LF, &
                               THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                               "test_read_expression_vectors_tsv_refuses_nan_and_inf: "//trim(WORDS(i_word)))
        end do

        tsv = scratch_file(scratch, "past_huge.tsv")
        call write_scratch_file(tsv, "id"//TAB//"a"//LF//"g1"//TAB//"1e999"//LF)
        file_list(1) = tsv
        call read_expression_vectors_tsv(file_list, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_true(is_err(ierr), "test_read_expression_vectors_tsv_refuses_nan_and_inf: 1e999")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_nan_and_inf

    !> Rows naming genes not in `gene_ids` are skipped however many there are: five rows for three
    !| genes, two of them unknown, are read. (More rows than genes used to be refused.)
    subroutine test_read_expression_vectors_tsv_row_count()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1)

        call create_scratch_directory(scratch, "expression_row_count")
        value_cols = [2]
        call read_one_file(scratch, "five_rows.tsv", "id"//TAB//"a"//LF//"gX"//TAB//"9"//LF// &
                           "g1"//TAB//"1"//LF//"g2"//TAB//"2"//LF//"gY"//TAB//"9"//LF//"g3"//TAB//"3"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_OK, &
                           "test_read_expression_vectors_tsv_row_count: five rows for three genes")
        expected = reshape([1.0_real64, 2.0_real64, 3.0_real64], [1, 3])
        call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_row_count: five rows for three genes")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_row_count

    !> Blank lines at the end of a file are no rows: an empty line, several, and one of blanks are
    !| read without an error. A blank line between two rows is a row without its fields and is
    !| refused with ERR_INVALID_INPUT naming file_list (argument 1).
    subroutine test_read_expression_vectors_tsv_trailing_empty_line()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3), expected(1, 3)
        integer(int32) :: value_cols(1), i_ending
        character(len=:), allocatable :: rows
        character(len=8), parameter :: ENDINGS(3) = ["one     ", "three   ", "blanks  "]

        call create_scratch_directory(scratch, "expression_trailing_empty_line")
        value_cols = [2]
        rows = "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF//"g2"//TAB//"2"//LF//"g3"//TAB//"3"//LF
        expected = reshape([1.0_real64, 2.0_real64, 3.0_real64], [1, 3])
        do i_ending = 1, size(ENDINGS)
            select case (i_ending)
            case (1)
                call read_one_file(scratch, "expression.tsv", rows//LF, THREE_GENE_IDS, value_cols, expression_vectors, &
                                   ERR_OK, "test_read_expression_vectors_tsv_trailing_empty_line: "//trim(ENDINGS(i_ending)))
            case (2)
                call read_one_file(scratch, "expression.tsv", rows//LF//LF//LF, THREE_GENE_IDS, value_cols, &
                                   expression_vectors, ERR_OK, &
                                   "test_read_expression_vectors_tsv_trailing_empty_line: "//trim(ENDINGS(i_ending)))
            case default
                call read_one_file(scratch, "expression.tsv", rows//"   "//LF, THREE_GENE_IDS, value_cols, &
                                   expression_vectors, ERR_OK, &
                                   "test_read_expression_vectors_tsv_trailing_empty_line: "//trim(ENDINGS(i_ending)))
            end select
            call assert_equal_array_real(expression_vectors, expected, 3, 0.0_real64, &
                                         "test_read_expression_vectors_tsv_trailing_empty_line: "//trim(ENDINGS(i_ending)))
        end do

        call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF//LF//"g3"//TAB//"3"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_trailing_empty_line: a blank line between rows")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_trailing_empty_line

    !> A gene named by two rows of one file is refused with ERR_INVALID_INPUT naming file_list
    !| (argument 1), with fewer rows than genes and with more. (The later row used to overwrite the
    !| earlier one.) In two files a gene has a row each.
    subroutine test_read_expression_vectors_tsv_refuses_duplicate_rows()
        type(scratch_directory_type) :: scratch
        real(real64) :: expression_vectors(1, 3)
        integer(int32) :: value_cols(1)

        call create_scratch_directory(scratch, "expression_duplicate_rows")
        value_cols = [2]
        call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF// &
                           "g1"//TAB//"1"//LF//"g2"//TAB//"2"//LF//"g1"//TAB//"7"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_refuses_duplicate_rows: three rows")
        call read_one_file(scratch, "expression.tsv", "id"//TAB//"a"//LF// &
                           "g1"//TAB//"1"//LF//"g2"//TAB//"2"//LF//"g3"//TAB//"3"//LF//"g1"//TAB//"4"//LF, &
                           THREE_GENE_IDS, value_cols, expression_vectors, ERR_INVALID_INPUT, &
                           "test_read_expression_vectors_tsv_refuses_duplicate_rows: four rows")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_duplicate_rows

    !> An empty cell is refused with ERR_INVALID_INPUT naming file_list (argument 1), wherever it
    !| is read: of "g1, 1, (empty), 3" the value column 3, and an empty id. An empty cell of a column
    !| that is not read is no error. (The cells used to be split with empty ones dropped, so the
    !| values right of an empty cell moved one column left.)
    subroutine test_read_expression_vectors_tsv_refuses_an_empty_cell()
        type(scratch_directory_type) :: scratch
        character(len=4) :: one_gene(1)
        real(real64) :: expression_vectors(2, 1), expected(2, 1)
        integer(int32) :: value_cols(2), ierr
        character(len=1024) :: file_list(1)

        call create_scratch_directory(scratch, "expression_empty_cell")
        file_list(1) = scratch_file(scratch, "expression.tsv")
        one_gene = ["g1"]
        value_cols = [2, 3]

        call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//TAB//"b"//TAB//"c"//LF//"g1"//TAB//"1"//TAB//TAB//"3"//LF)
        call read_expression_vectors_tsv(file_list, one_gene, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_refuses_an_empty_cell: a value", 1)

        call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//TAB//"b"//LF//TAB//"1"//TAB//"2"//LF)
        call read_expression_vectors_tsv(file_list, one_gene, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_refuses_an_empty_cell: the id", 1)

        call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//TAB//"b"//TAB//"c"//LF//"g1"//TAB//"1"//TAB//"2"//TAB//LF)
        expression_vectors = UNWRITTEN
        call read_expression_vectors_tsv(file_list, one_gene, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_refuses_an_empty_cell: an unread one, ierr")
        expected = reshape([1.0_real64, 2.0_real64], [2, 1])
        call assert_equal_array_real(expression_vectors, expected, 2, 0.0_real64, &
                                     "test_read_expression_vectors_tsv_refuses_an_empty_cell: an unread one")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_an_empty_cell

    !> A cell that is not one complete number is refused with ERR_READ_DATA naming file_list
    !| (argument 1): a decimal comma (0,25), a slash, two numbers (1 2, 1;2) and a repeat count
    !| (2*3). They used to be read list-directed, which stops at a comma or a blank and takes a slash
    !| for "no value": 0,25 read as 0, and a slash left the previous value.
    subroutine test_read_expression_vectors_tsv_refuses_a_partial_number()
        type(scratch_directory_type) :: scratch
        character(len=4) :: one_gene(1)
        real(real64) :: expression_vectors(2, 1)
        integer(int32) :: value_cols(2), ierr, i_cell
        character(len=1024) :: file_list(1)
        character(len=4), parameter :: CELLS(5) = ["0,25", "/   ", "1 2 ", "1;2 ", "2*3 "]

        call create_scratch_directory(scratch, "expression_partial_number")
        file_list(1) = scratch_file(scratch, "expression.tsv")
        one_gene = ["g1"]
        value_cols = [2, 3]
        do i_cell = 1, size(CELLS)
            call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//TAB//"b"//LF//"g1"//TAB//"1.5"//TAB// &
                                    trim(CELLS(i_cell))//LF)
            call read_expression_vectors_tsv(file_list, one_gene, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
            call assert_err(ierr, ERR_READ_DATA, "test_read_expression_vectors_tsv_refuses_a_partial_number: " &
                            //trim(CELLS(i_cell)), 1)
        end do

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_refuses_a_partial_number

    !> A row's id is looked up whole: with ids of 8 characters, a row of NP_001.10 names no gene of
    !| the list and is skipped, without an error, and NP_001.1's column is left alone. (The id used
    !| to be cut to 8 characters first, and the row read into NP_001.1's column.)
    subroutine test_read_expression_vectors_tsv_matches_whole_ids()
        type(scratch_directory_type) :: scratch
        character(len=8) :: gene_ids(1)
        real(real64) :: expression_vectors(1, 1)
        integer(int32) :: value_cols(1), ierr
        character(len=1024) :: file_list(1)

        call create_scratch_directory(scratch, "expression_whole_ids")
        file_list(1) = scratch_file(scratch, "expression.tsv")
        call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//LF//"NP_001.10"//TAB//"5"//LF)
        gene_ids = ["NP_001.1"]
        value_cols = [2]
        expression_vectors = UNWRITTEN

        call read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_matches_whole_ids: ierr")
        call assert_equal_real(expression_vectors(1, 1), UNWRITTEN, 0.0_real64, &
                               "test_read_expression_vectors_tsv_matches_whole_ids: NP_001.1's column")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_matches_whole_ids

    !> A line is read whole up to 4096 characters (the readers' limit) and refused with
    !| ERR_STRING_TOO_LONG naming file_list (argument 1) past it: 9.75 at characters 2047-2050
    !| (past the 2048 characters this reader used to cut lines at, silently), and at 4093-4096, is
    !| read; one character more is refused.
    subroutine test_read_expression_vectors_tsv_reads_a_long_line_whole()
        type(scratch_directory_type) :: scratch
        character(len=4096) :: gene_ids(1)
        real(real64) :: expression_vectors(1, 1)
        integer(int32) :: value_cols(1), ierr, i_length
        character(len=1024) :: file_list(1)
        integer(int32), parameter :: ID_LENGTHS(3) = [2045, 4091, 4092]

        call create_scratch_directory(scratch, "expression_long_line")
        file_list(1) = scratch_file(scratch, "expression.tsv")
        value_cols = [2]
        do i_length = 1, size(ID_LENGTHS)
            ! the id, a tab, and 9.75: a line of the id's length + 5
            gene_ids(1) = repeat("x", ID_LENGTHS(i_length))
            call write_scratch_file(trim(file_list(1)), "id"//TAB//"a"//LF//trim(gene_ids(1))//TAB//"9.75"//LF)
            expression_vectors = UNWRITTEN
            call read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
            if (ID_LENGTHS(i_length) + 5 <= 4096) then
                call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_reads_a_long_line_whole: line of " &
                                //(ID_LENGTHS(i_length) + 5)//", ierr")
                call assert_equal_real(expression_vectors(1, 1), 9.75_real64, 0.0_real64, &
                                       "test_read_expression_vectors_tsv_reads_a_long_line_whole: line of " &
                                       //(ID_LENGTHS(i_length) + 5))
            else
                call assert_err(ierr, ERR_STRING_TOO_LONG, "test_read_expression_vectors_tsv_reads_a_long_line_whole: line of " &
                                //(ID_LENGTHS(i_length) + 5), 1)
            end if
        end do

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_reads_a_long_line_whole

    !> A file that does not exist is refused with ERR_FILE_OPEN, as the only file and as the
    !| second after one that was read, and the error unit stays connected (see
    !| `test_read_gene_ids_from_tsv_file_missing_file`).
    subroutine test_read_expression_vectors_tsv_missing_file()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: one_file(1), two_files(2)
        real(real64) :: expression_vectors(2, 3)
        integer(int32) :: value_cols(1), ierr

        call create_scratch_directory(scratch, "expression_missing_file")
        one_file(1) = scratch_file(scratch, "missing.tsv")
        two_files(1) = scratch_file(scratch, "present.tsv")
        call write_scratch_file(trim(two_files(1)), "id"//TAB//"a"//LF//"g1"//TAB//"1"//LF)
        two_files(2) = one_file(1)
        value_cols = [2]

        call clear_stack()
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "test_read_expression_vectors_tsv_missing_file: the only file")
        call check_error_unit_connected("test_read_expression_vectors_tsv_missing_file: the only file")
        call clear_stack()
        call read_expression_vectors_tsv(two_files, THREE_GENE_IDS, expression_vectors, 1, 1, value_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "test_read_expression_vectors_tsv_missing_file: the second file")
        call check_error_unit_connected("test_read_expression_vectors_tsv_missing_file: the second file")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_missing_file

    !> `gene_col` (argument 5), every entry of `value_cols` (argument 6) and `start_row`
    !| (argument 7) count from 1: 0 is refused with ERR_INVALID_INPUT naming the argument, and so is
    !| a -1 among value_cols. And `expression_vectors` (argument 3) needs a row for every value
    !| column of every file from start_row on: two columns into one row, start_row 2 for two
    !| columns in two rows, and two files of one column into one row are refused with
    !| ERR_SIZE_MISMATCH; start_row 1 for two columns in two rows is read. They used to read or
    !| write past the arrays' ends.
    subroutine test_read_expression_vectors_tsv_checks_the_positions()
        type(scratch_directory_type) :: scratch
        character(len=1024) :: one_file(1), two_files(2)
        real(real64) :: one_row(1, 3), two_rows(2, 3)
        integer(int32) :: two_cols(2), one_col(1), with_zero(2), with_negative(2), ierr

        call create_scratch_directory(scratch, "expression_positions")
        one_file(1) = scratch_file(scratch, "expression.tsv")
        call write_scratch_file(trim(one_file(1)), "id"//TAB//"a"//TAB//"b"//LF//"g1"//TAB//"1"//TAB//"2"//LF)
        two_files = one_file(1)
        two_cols = [2, 3]
        one_col = [2]
        with_zero = [2, 0]
        with_negative = [-1, 2]

        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 0, two_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_checks_the_positions: gene_col 0", 5)
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 1, with_zero, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_checks_the_positions: value_cols 0", 6)
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 1, with_negative, 1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_checks_the_positions: value_cols -1", 6)
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 1, two_cols, 0, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_read_expression_vectors_tsv_checks_the_positions: start_row 0", 7)

        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, one_row, 1, 1, two_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_expression_vectors_tsv_checks_the_positions: 2 columns, 1 row", 3)
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 1, two_cols, 2, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, &
                        "test_read_expression_vectors_tsv_checks_the_positions: start_row 2, 2 columns, 2 rows", 3)
        call read_expression_vectors_tsv(two_files, THREE_GENE_IDS, one_row, 1, 1, one_col, 1, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_expression_vectors_tsv_checks_the_positions: 2 files, 1 row", 3)
        call read_expression_vectors_tsv(one_file, THREE_GENE_IDS, two_rows, 1, 1, two_cols, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_read_expression_vectors_tsv_checks_the_positions: start_row 1, 2 rows")

        call remove_scratch_directory(scratch)
    end subroutine test_read_expression_vectors_tsv_checks_the_positions
end module mod_test_read_expression_vectors_tsv
