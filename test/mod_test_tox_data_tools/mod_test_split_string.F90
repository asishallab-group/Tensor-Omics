!> The `split_string` cases: the fields of a line between delimiters, each without its surrounding
!| blanks; the default delimiter, a blank, and degenerate lines (empty, nothing but delimiters,
!| one field). The length of the strings that hold the fields is not documented and not pinned.
!|
!| How an empty field between two delimiters other than a blank is split is not pinned here: the
!| TSV reader's handling of an empty cell rests on it, and that is an open bug
!| (`test_read_expression_vectors_tsv_refuses_an_empty_cell`).
module mod_test_split_string
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_tools, only: split_string
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_data_tools_fixtures, only: TAB
    implicit none
    private
    public :: get_all_tests_split_string

contains

    !> Every case of `split_string`.
    function get_all_tests_split_string() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(5))
        all_tests(1) = test_case("test_split_string_fields", test_split_string_fields)
        all_tests(2) = test_case("test_split_string_trims_the_fields", test_split_string_trims_the_fields)
        all_tests(3) = test_case("test_split_string_default_delimiter", test_split_string_default_delimiter)
        all_tests(4) = test_case("test_split_string_one_field", test_split_string_one_field)
        all_tests(5) = test_case("test_split_string_no_field", test_split_string_no_field)
    end function get_all_tests_split_string

    !> Splits `input` at `delimiter` and checks ierr OK, the number of fields and each field.
    subroutine check_split(input, delimiter, expected, msg)
        character(*), intent(in) :: input, expected(:), msg
        character(len=1), intent(in) :: delimiter
        character(len=:), allocatable :: fields(:)
        integer(int32) :: ierr, i_field

        call split_string(input, fields, delimiter, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_int(size(fields), size(expected), msg//": number of fields")
        if (size(fields) /= size(expected)) return
        do i_field = 1, size(expected)
            call assert_string_equal(fields(i_field), expected(i_field), msg//": field "//i_field)
        end do
    end subroutine check_split

    !> Three tab-separated fields of 1, 2 and 3 characters; a comma is part of a field when the
    !| delimiter is a tab, and the delimiter when it is a comma.
    subroutine test_split_string_fields()
        character(len=4) :: expected(3)

        expected = ["a  ", "bb ", "ccc"]
        call check_split("a"//TAB//"bb"//TAB//"ccc", TAB, expected, "test_split_string_fields: tab")
        expected = ["x,y", "z  ", "1  "]
        call check_split("x,y"//TAB//"z"//TAB//"1", TAB, expected, "test_split_string_fields: a comma in a tab field")
        expected = ["x", "y", "z"]
        call check_split("x,y,z", ",", expected, "test_split_string_fields: comma")
    end subroutine test_split_string_fields

    !> Blanks around a field are not part of it, so "a, b ,c" splits at commas into a, b and c.
    subroutine test_split_string_trims_the_fields()
        character(len=1) :: expected(3)

        expected = ["a", "b", "c"]
        call check_split("a, b ,c", ",", expected, "test_split_string_trims_the_fields")
    end subroutine test_split_string_trims_the_fields

    !> Without a delimiter the fields are separated by blanks, and a run of blanks, leading and
    !| trailing ones included, separates once: "  ab  c d  " is ab, c and d. A tab is no blank.
    subroutine test_split_string_default_delimiter()
        character(len=:), allocatable :: fields(:)
        integer(int32) :: ierr

        call split_string("  ab  c d  ", fields, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_split_string_default_delimiter: ierr")
        call assert_equal_int(size(fields), 3, "test_split_string_default_delimiter: number of fields")
        if (size(fields) == 3) then
            call assert_string_equal(fields(1), "ab", "test_split_string_default_delimiter: field 1")
            call assert_string_equal(fields(2), "c", "test_split_string_default_delimiter: field 2")
            call assert_string_equal(fields(3), "d", "test_split_string_default_delimiter: field 3")
        end if

        call split_string("a"//TAB//"b", fields, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_split_string_default_delimiter: a tab, ierr")
        call assert_equal_int(size(fields), 1, "test_split_string_default_delimiter: a tab is no blank")
    end subroutine test_split_string_default_delimiter

    !> A line without the delimiter is one field, the line.
    subroutine test_split_string_one_field()
        character(len=5) :: expected(1)

        expected = ["gene1"]
        call check_split("gene1", TAB, expected, "test_split_string_one_field")
    end subroutine test_split_string_one_field

    !> An empty line, a blank one, and one of nothing but delimiters have no field: an array of
    !| none, without an error.
    subroutine test_split_string_no_field()
        character(len=1) :: expected(0)

        call check_split("", TAB, expected, "test_split_string_no_field: empty")
        call check_split("   ", ",", expected, "test_split_string_no_field: blank")
        call check_split(",,,", ",", expected, "test_split_string_no_field: only delimiters")
    end subroutine test_split_string_no_field
end module mod_test_split_string
