!> Unit tests for the streaming JSON writer, f42_serde_json_serialize.
!|
!| Every case writes its own file into the working directory and removes it again. A file is
!| compared byte for byte only where every compiler must produce the same bytes (integers,
!| strings, and reals whose 17-digit form is exact); other reals are read back and compared bit
!| for bit, because compilers may round the 17th digit differently.
!|
!| The file helpers at the end are public: the flyer suite uses them too.
module mod_test_json_serialize
    use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use asserts
    use test_suite
    use tox_errors, only: ERR_OK, ERR_FILE_OPEN, ERR_JSON_SEQUENCE, ERR_NAN_INF, ERR_INVALID_UTF8, ERR_INVALID_INPUT
    use f42_serde_json_serialize
    implicit none
    private
    public :: get_all_tests_json_serialize
    public :: read_whole_file, file_exists, remove_file, assert_file_equals, assert_no_file

    character(len=*), parameter :: LF = achar(10)
    character(len=*), parameter :: BS = achar(92)

contains

    function get_all_tests_json_serialize() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(40))
        all_tests(1) = test_case("test_empty_object_root", test_empty_object_root)
        all_tests(2) = test_case("test_empty_array_root", test_empty_array_root)
        all_tests(3) = test_case("test_scalar_roots", test_scalar_roots)
        all_tests(4) = test_case("test_nested_document", test_nested_document)
        all_tests(5) = test_case("test_integer_extremes", test_integer_extremes)
        all_tests(6) = test_case("test_integer_arrays_across_blocks", test_integer_arrays_across_blocks)
        all_tests(7) = test_case("test_real_exact_text", test_real_exact_text)
        all_tests(8) = test_case("test_real_round_trip", test_real_round_trip)
        all_tests(9) = test_case("test_real_array_across_blocks", test_real_array_across_blocks)
        all_tests(10) = test_case("test_nan_refused", test_nan_refused)
        all_tests(11) = test_case("test_infinity_refused_in_array", test_infinity_refused_in_array)
        all_tests(12) = test_case("test_string_escapes", test_string_escapes)
        all_tests(13) = test_case("test_string_trailing_blanks", test_string_trailing_blanks)
        all_tests(14) = test_case("test_utf8_valid_passthrough", test_utf8_valid_passthrough)
        all_tests(15) = test_case("test_utf8_invalid_refused", test_utf8_invalid_refused)
        all_tests(16) = test_case("test_utf8_validator", test_utf8_validator)
        all_tests(17) = test_case("test_is_finite", test_is_finite)
        all_tests(18) = test_case("test_value_without_key_in_object", test_value_without_key_in_object)
        all_tests(19) = test_case("test_key_in_array", test_key_in_array)
        all_tests(20) = test_case("test_key_at_root", test_key_at_root)
        all_tests(21) = test_case("test_mismatched_end", test_mismatched_end)
        all_tests(22) = test_case("test_end_without_container", test_end_without_container)
        all_tests(23) = test_case("test_second_root", test_second_root)
        all_tests(24) = test_case("test_unclosed_container_at_close", test_unclosed_container_at_close)
        all_tests(25) = test_case("test_no_root_at_close", test_no_root_at_close)
        all_tests(26) = test_case("test_errors_are_sticky", test_errors_are_sticky)
        all_tests(27) = test_case("test_writer_not_open", test_writer_not_open)
        all_tests(28) = test_case("test_open_twice", test_open_twice)
        all_tests(29) = test_case("test_writer_reuse", test_writer_reuse)
        all_tests(30) = test_case("test_existing_file_refused", test_existing_file_refused)
        all_tests(31) = test_case("test_unwritable_path", test_unwritable_path)
        all_tests(32) = test_case("test_deep_nesting", test_deep_nesting)
        all_tests(33) = test_case("test_string_longer_than_buffer", test_string_longer_than_buffer)
        all_tests(34) = test_case("test_output_larger_than_buffer", test_output_larger_than_buffer)
        all_tests(35) = test_case("test_fuzz_call_sequences", test_fuzz_call_sequences)
        all_tests(36) = test_case("test_boolean_kinds", test_boolean_kinds)
        all_tests(37) = test_case("test_unsupported_types", test_unsupported_types)
        all_tests(38) = test_case("test_zero_decided_from_bits", test_zero_decided_from_bits)
        all_tests(39) = test_case("test_empty_arrays", test_empty_arrays)
        all_tests(40) = test_case("test_array_sections", test_array_sections)
    end function get_all_tests_json_serialize

    ! ============================================================================================
    ! Documents
    ! ============================================================================================

    subroutine test_empty_object_root()
        character(len=*), parameter :: filename = "json_empty_object.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call assert_err(ierr, ERR_OK, "open succeeds")
        call json_begin_object(w)
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "an empty object is a document")
        call assert_file_equals(filename, "{}"//LF, "empty object")
        call remove_file(filename)
    end subroutine test_empty_object_root

    subroutine test_empty_array_root()
        character(len=*), parameter :: filename = "json_empty_array.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "an empty array is a document")
        call assert_file_equals(filename, "[]"//LF, "empty array")
        call remove_file(filename)
    end subroutine test_empty_array_root

    !> Any value may be the root.
    subroutine test_scalar_roots()
        character(len=*), parameter :: filename = "json_scalar_root.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, i_case

        do i_case = 1, 6
            call remove_file(filename)
            call json_open(w, filename, ierr)
            select case (i_case)
            case (1)
                call json_element(w, 42_int32)
            case (2)
                call json_element(w, -7_int64)
            case (3)
                call json_element(w, 0.5_real64)
            case (4)
                call json_element(w, "root")
            case (5)
                call json_element(w, .false.)
            case (6)
                call json_null(w)
            end select
            call json_close(w, ierr)
            call assert_err(ierr, ERR_OK, "scalar root "//i_case)
            select case (i_case)
            case (1)
                call assert_file_equals(filename, "42"//LF, "int32 root")
            case (2)
                call assert_file_equals(filename, "-7"//LF, "int64 root")
            case (3)
                call assert_file_equals(filename, "5.0000000000000000E-001"//LF, "real root")
            case (4)
                call assert_file_equals(filename, '"root"'//LF, "string root")
            case (5)
                call assert_file_equals(filename, "false"//LF, "logical root")
            case (6)
                call assert_file_equals(filename, "null"//LF, "null root")
            end select
        end do
        call remove_file(filename)
    end subroutine test_scalar_roots

    !> Every construct once: members of every type, arrays of every type, empty and nested containers.
    subroutine test_nested_document()
        character(len=*), parameter :: filename = "json_nested.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr
        logical(c_bool), parameter :: flags(3) = [.true._c_bool, .false._c_bool, .true._c_bool]
        character(len=5), parameter :: words(3) = ["ab   ", "c    ", "     "]
        integer(int32) :: no_values(0)
        character(len=:), allocatable :: expected

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, "i", 1_int32)
        call json_member(w, "l", 2_int64)
        call json_member(w, "r", -2.0_real64)
        call json_member(w, "b", .true._c_bool)
        call json_member(w, "d", .false.)
        call json_member(w, "s", "text")
        call json_null(w, "n")
        call json_member(w, "ia", [1_int32, -2_int32, 3_int32])
        call json_member(w, "la", [4_int64])
        call json_member(w, "ra", [0.25_real64, 0.0_real64])
        call json_member(w, "ba", flags)
        call json_member(w, "da", [.false., .true.])
        call json_member(w, "sa", words)
        call json_member(w, "empty", no_values)
        call json_begin_object(w, "o")
        call json_end_object(w)
        call json_begin_array(w, "a")
        call json_element(w, 5_int32)
        call json_element(w, [6_int64, 7_int64])
        call json_begin_object(w)
        call json_member(w, "x", "y")
        call json_end_object(w)
        call json_begin_array(w)
        call json_end_array(w)
        call json_null(w)
        call json_element(w, .true._c_bool)
        call json_element(w, 8.0_real64)
        call json_element(w, "z")
        call json_end_array(w)
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "nested document")

        expected = '{"i":1,"l":2,"r":-2.0000000000000000E+000,"b":true,"d":false,"s":"text","n":null,' &
                   //'"ia":[1,-2,3],"la":[4],"ra":[2.5000000000000000E-001,0.0],"ba":[true,false,true],' &
                   //'"da":[false,true],"sa":["ab","c",""],"empty":[],"o":{},' &
                   //'"a":[5,[6,7],{"x":"y"},[],null,true,8.0000000000000000E+000,"z"]}'//LF
        call assert_file_equals(filename, expected, "nested document")
        call remove_file(filename)
    end subroutine test_nested_document

    ! ============================================================================================
    ! Numbers
    ! ============================================================================================

    subroutine test_integer_extremes()
        character(len=*), parameter :: filename = "json_integers.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, lowest_int32, pair32(2)
        integer(int64) :: lowest_int64, pair64(2)

        ! the most negative values, computed because no literal can spell them
        lowest_int32 = -huge(1_int32)
        lowest_int32 = lowest_int32 - 1_int32
        lowest_int64 = -huge(1_int64)
        lowest_int64 = lowest_int64 - 1_int64
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, 0_int32)
        call json_element(w, -1_int32)
        call json_element(w, huge(1_int32))
        call json_element(w, lowest_int32)
        call json_element(w, huge(1_int64))
        call json_element(w, lowest_int64)
        pair32(1) = huge(1_int32)
        pair32(2) = lowest_int32
        pair64(1) = huge(1_int64)
        pair64(2) = lowest_int64
        call json_element(w, pair32)
        call json_element(w, pair64)
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "integer extremes")
        call assert_file_equals(filename, "[0,-1,2147483647,-2147483648,9223372036854775807,-9223372036854775808," &
                                //"[2147483647,-2147483648],[9223372036854775807,-9223372036854775808]]"//LF, &
                                "integer extremes")
        call remove_file(filename)
    end subroutine test_integer_extremes

    !> Arrays longer than one formatting block keep every element and one comma between neighbours.
    subroutine test_integer_arrays_across_blocks()
        character(len=*), parameter :: filename = "json_integer_blocks.test.json"
        integer(int32), parameter :: n_values = 600
        type(json_writer) :: w
        integer(int32) :: ierr, i_value
        integer(int32) :: values32(n_values)
        integer(int64) :: values64(n_values)
        character(len=:), allocatable :: expected32, expected64
        character(len=24) :: digits

        do i_value = 1, n_values
            values32(i_value) = i_value*(-1)**i_value
            values64(i_value) = int(i_value, int64)*10000000000_int64
        end do
        expected32 = "["
        expected64 = "["
        do i_value = 1, n_values
            if (i_value > 1) then
                expected32 = expected32//","
                expected64 = expected64//","
            end if
            write (digits, "(I0)") values32(i_value)
            expected32 = expected32//trim(digits)
            write (digits, "(I0)") values64(i_value)
            expected64 = expected64//trim(digits)
        end do

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, values32)
        call json_element(w, values64)
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "integer arrays")
        call assert_file_equals(filename, "["//expected32//"],"//expected64//"]]"//LF, "integer arrays across blocks")
        call remove_file(filename)
    end subroutine test_integer_arrays_across_blocks

    !> Reals whose 17 significant digits are exact print identically everywhere; zeros come from the bits.
    subroutine test_real_exact_text()
        character(len=*), parameter :: filename = "json_real_text.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, 0.0_real64)
        call json_element(w, -0.0_real64)
        call json_element(w, 1.0_real64)
        call json_element(w, -2.0_real64)
        call json_element(w, 0.5_real64)
        call json_element(w, 1024.0_real64)
        call json_element(w, -0.125_real64)
        call json_element(w, 2.0_real64**(-20))
        call json_element(w, [0.0_real64, -0.0_real64, 3.0_real64])
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "exact reals")
        call assert_file_equals(filename, "[0.0,-0.0,1.0000000000000000E+000,-2.0000000000000000E+000," &
                                //"5.0000000000000000E-001,1.0240000000000000E+003,-1.2500000000000000E-001," &
                                //"9.5367431640625000E-007,[0.0,-0.0,3.0000000000000000E+000]]"//LF, "exact reals")
        call remove_file(filename)
    end subroutine test_real_exact_text

    !> 17 digits read back as the same double, including the extremes and subnormals.
    subroutine test_real_round_trip()
        character(len=*), parameter :: filename = "json_real_round_trip.test.json"
        integer(int32), parameter :: n_values = 14
        type(json_writer) :: w
        integer(int32) :: ierr, i_value
        real(real64) :: values(n_values), parsed(n_values)

        values = [0.1_real64, 1.0_real64/3.0_real64, acos(-1.0_real64), huge(1.0_real64), -huge(1.0_real64), &
                  tiny(1.0_real64), -tiny(1.0_real64), transfer(1_int64, 1.0_real64), &
                  transfer(int(z'000FFFFFFFFFFFFF', int64), 1.0_real64), 1.0e300_real64, -1.0e-300_real64, &
                  -0.0_real64, 123456789.123456789_real64, nearest(1.0_real64, 2.0_real64)]

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_element(w, values)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "round trip written")
        call read_real_array_file(filename, parsed, ierr)
        call assert_equal_int(ierr, 0, "round trip file parses")
        do i_value = 1, n_values
            call assert_true(transfer(parsed(i_value), 0_int64) == transfer(values(i_value), 0_int64), &
                             "value "//i_value//" reads back bit for bit")
        end do
        call remove_file(filename)
    end subroutine test_real_round_trip

    !> An array spanning three formatting blocks, zeros at the block edges included.
    subroutine test_real_array_across_blocks()
        character(len=*), parameter :: filename = "json_real_blocks.test.json"
        integer(int32), parameter :: n_values = 600
        type(json_writer) :: w
        integer(int32) :: ierr, i_value
        real(real64) :: values(n_values), parsed(n_values)

        do i_value = 1, n_values
            values(i_value) = sqrt(real(i_value, real64))*(-1.0_real64)**i_value
        end do
        values([1, 256, 257, 512, 513, 600]) = [0.0_real64, -0.0_real64, 0.0_real64, 0.0_real64, -0.0_real64, 0.0_real64]

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_element(w, values)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "real blocks written")
        call read_real_array_file(filename, parsed, ierr)
        call assert_equal_int(ierr, 0, "real blocks parse")
        call assert_equal_int(count_bit_mismatches(parsed, values), 0, "every value of three blocks reads back bit for bit")
        call remove_file(filename)
    end subroutine test_real_array_across_blocks

    subroutine test_nan_refused()
        character(len=*), parameter :: filename = "json_nan.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, "x", ieee_value(1.0_real64, ieee_quiet_nan))
        call assert_err(w%ierr, ERR_NAN_INF, "NaN member")
        call assert_equal_int(int(w%failing_call, int32), 3, "NaN fails the third call")
        call assert_no_file(filename, "NaN deletes the file at once")
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_NAN_INF, "close returns the NaN error")
        call assert_no_file(filename, "NaN leaves no file")
    end subroutine test_nan_refused

    subroutine test_infinity_refused_in_array()
        character(len=*), parameter :: filename = "json_infinity.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, i_case
        real(real64) :: values(300)

        do i_case = 1, 3
            values = 1.0_real64
            select case (i_case)
            case (1)
                values(300) = ieee_value(1.0_real64, ieee_positive_inf)
            case (2)
                values(1) = ieee_value(1.0_real64, ieee_negative_inf)
            case (3)
                values(257) = ieee_value(1.0_real64, ieee_quiet_nan)
            end select
            call remove_file(filename)
            call json_open(w, filename, ierr)
            call json_begin_array(w)
            call json_element(w, values)
            call json_close(w, ierr)
            call assert_err(ierr, ERR_NAN_INF, "non-finite array element, case "//i_case)
            call assert_equal_int(int(w%failing_call, int32), 3, "the array call fails, case "//i_case)
            call assert_no_file(filename, "non-finite element leaves no file, case "//i_case)
        end do
    end subroutine test_infinity_refused_in_array

    ! ============================================================================================
    ! Strings
    ! ============================================================================================

    !> Every byte below 0x20, the quote and the backslash are escaped; / and DEL are not.
    subroutine test_string_escapes()
        character(len=*), parameter :: filename = "json_escapes.test.json"
        character(len=*), parameter :: HEX = "0123456789abcdef"
        type(json_writer) :: w
        integer(int32) :: ierr, code
        character(len=35) :: text
        character(len=:), allocatable :: expected

        do code = 0, 31
            text(code + 1:code + 1) = achar(code)
        end do
        text(33:35) = '"'//BS//"/"
        expected = '"'
        do code = 0, 31
            select case (code)
            case (8)
                expected = expected//BS//"b"
            case (9)
                expected = expected//BS//"t"
            case (10)
                expected = expected//BS//"n"
            case (12)
                expected = expected//BS//"f"
            case (13)
                expected = expected//BS//"r"
            case default
                expected = expected//BS//"u00"//HEX(code/16 + 1:code/16 + 1)//HEX(mod(code, 16) + 1:mod(code, 16) + 1)
            end select
        end do
        expected = expected//BS//'"'//BS//BS//"/"//char(127)//'"'

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, text//char(127), text//char(127))
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "escapes")
        call assert_file_equals(filename, "{"//expected//":"//expected//"}"//LF, "keys and values escaped alike")
        call remove_file(filename)
    end subroutine test_string_escapes

    subroutine test_string_trailing_blanks()
        character(len=*), parameter :: filename = "json_blanks.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, "key   ", "abc   ")
        call json_member(w, "  lead", "  lead  ")
        call json_member(w, "blank", "     ")
        call json_member(w, "", "")
        call json_member(w, "tab", "a"//achar(9)//"  ")
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "blanks")
        call assert_file_equals(filename, '{"key":"abc","  lead":"  lead","blank":"","":"","tab":"a'//BS//'t"}'//LF, &
                                "only trailing blanks are removed")
        call remove_file(filename)
    end subroutine test_string_trailing_blanks

    subroutine test_utf8_valid_passthrough()
        character(len=*), parameter :: filename = "json_utf8_valid.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr
        character(len=:), allocatable :: text

        ! U+0080, U+00E9, U+07FF, U+0800, U+20AC, U+D7FF, U+E000, U+FFFF, U+10000, U+1F600, U+10FFFF
        text = bytes([194, 128]) // bytes([195, 169]) // bytes([223, 191]) // bytes([224, 160, 128]) &
               //bytes([226, 130, 172]) // bytes([237, 159, 191]) // bytes([238, 128, 128]) &
               //bytes([239, 191, 191]) // bytes([240, 144, 128, 128]) // bytes([240, 159, 152, 128]) &
               //bytes([244, 143, 191, 191])

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, text, text//"x")
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "valid UTF-8")
        call assert_file_equals(filename, '{"'//text//'":"'//text//'x"}'//LF, "valid UTF-8 copied unchanged")
        call remove_file(filename)
    end subroutine test_utf8_valid_passthrough

    !> Each malformed sequence on its own, as a value and as a key.
    subroutine test_utf8_invalid_refused()
        character(len=*), parameter :: filename = "json_utf8_invalid.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, i_case, i_role
        character(len=8) :: bad

        do i_case = 1, n_invalid_utf8_cases()
            bad = invalid_utf8_case(i_case)
            do i_role = 1, 2
                call remove_file(filename)
                call json_open(w, filename, ierr)
                call json_begin_object(w)
                if (i_role == 1) then
                    call json_member(w, "k", "ok"//trim(bad))
                else
                    call json_member(w, "ok"//trim(bad), "v")
                end if
                call json_end_object(w)
                call json_close(w, ierr)
                call assert_err(ierr, ERR_INVALID_UTF8, "invalid UTF-8 case "//i_case//" role "//i_role)
                call assert_equal_int(int(w%failing_call, int32), 3, "member call fails, case "//i_case)
                call assert_no_file(filename, "invalid UTF-8 leaves no file, case "//i_case)
            end do
        end do
    end subroutine test_utf8_invalid_refused

    subroutine test_utf8_validator()
        integer(int32) :: i_case

        call assert_true(json_is_valid_utf8(""), "empty is valid")
        call assert_true(json_is_valid_utf8("plain ascii "//achar(0)//char(127)), "ASCII is valid")
        call assert_true(json_is_valid_utf8(bytes([240, 159, 152, 128])//"   "), "emoji then blanks")
        call assert_true(json_is_valid_utf8(bytes([244, 143, 191, 191])), "U+10FFFF")
        do i_case = 1, n_invalid_utf8_cases()
            call assert_false(json_is_valid_utf8("a"//invalid_utf8_case(i_case)), "invalid case "//i_case)
        end do
    end subroutine test_utf8_validator

    subroutine test_is_finite()
        call assert_false(json_is_finite(ieee_value(1.0_real64, ieee_quiet_nan)), "NaN")
        call assert_false(json_is_finite(ieee_value(1.0_real64, ieee_positive_inf)), "+Inf")
        call assert_false(json_is_finite(ieee_value(1.0_real64, ieee_negative_inf)), "-Inf")
        call assert_true(json_is_finite(huge(1.0_real64)), "huge")
        call assert_true(json_is_finite(-huge(1.0_real64)), "-huge")
        call assert_true(json_is_finite(transfer(1_int64, 1.0_real64)), "smallest subnormal")
        call assert_true(json_is_finite(-0.0_real64), "-0")
    end subroutine test_is_finite

    ! ============================================================================================
    ! The state machine: each violation, its call number, and no file left
    ! ============================================================================================

    subroutine test_value_without_key_in_object()
        character(len=*), parameter :: filename = "json_seq_nokey.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_begin_object(w)
        call json_element(w, 1_int32)
        call assert_sequence_failure(w, filename, 3, "value without key in an object")

        call open_fresh(w, filename)
        call json_begin_object(w)
        call json_begin_array(w)
        call assert_sequence_failure(w, filename, 3, "array without key in an object")

        call open_fresh(w, filename)
        call json_begin_object(w)
        call json_null(w)
        call assert_sequence_failure(w, filename, 3, "null without key in an object")
    end subroutine test_value_without_key_in_object

    subroutine test_key_in_array()
        character(len=*), parameter :: filename = "json_seq_keyarray.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_member(w, "k", 1_int32)
        call assert_sequence_failure(w, filename, 3, "member in an array")

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_element(w, 1_int32)
        call json_null(w, "k")
        call assert_sequence_failure(w, filename, 4, "keyed null in an array")

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_begin_object(w, "k")
        call assert_sequence_failure(w, filename, 3, "keyed object in an array")
    end subroutine test_key_in_array

    subroutine test_key_at_root()
        character(len=*), parameter :: filename = "json_seq_keyroot.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_member(w, "k", "v")
        call assert_sequence_failure(w, filename, 2, "member as the root")

        call open_fresh(w, filename)
        call json_begin_array(w, "k")
        call assert_sequence_failure(w, filename, 2, "keyed array as the root")
    end subroutine test_key_at_root

    subroutine test_mismatched_end()
        character(len=*), parameter :: filename = "json_seq_mismatch.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_begin_object(w)
        call json_end_array(w)
        call assert_sequence_failure(w, filename, 3, "object ended as an array")

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_begin_array(w)
        call json_end_array(w)
        call json_end_object(w)
        call assert_sequence_failure(w, filename, 5, "array ended as an object")
    end subroutine test_mismatched_end

    subroutine test_end_without_container()
        character(len=*), parameter :: filename = "json_seq_noopen.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_end_object(w)
        call assert_sequence_failure(w, filename, 2, "end with nothing open")

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_end_array(w)
        call json_end_array(w)
        call assert_sequence_failure(w, filename, 4, "one end too many")
    end subroutine test_end_without_container

    subroutine test_second_root()
        character(len=*), parameter :: filename = "json_seq_tworoots.test.json"
        type(json_writer) :: w

        call open_fresh(w, filename)
        call json_element(w, 1_int32)
        call json_element(w, 2_int32)
        call assert_sequence_failure(w, filename, 3, "second scalar root")

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_end_array(w)
        call json_begin_object(w)
        call assert_sequence_failure(w, filename, 4, "second container root")
    end subroutine test_second_root

    subroutine test_unclosed_container_at_close()
        character(len=*), parameter :: filename = "json_seq_unclosed.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call open_fresh(w, filename)
        call json_begin_object(w)
        call json_begin_array(w, "a")
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "close with an open object")
        call assert_equal_int(int(w%failing_call, int32), 5, "close is the failing call")
        call assert_no_file(filename, "unclosed object leaves no file")
    end subroutine test_unclosed_container_at_close

    subroutine test_no_root_at_close()
        character(len=*), parameter :: filename = "json_seq_noroot.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call open_fresh(w, filename)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "close without a root")
        call assert_equal_int(int(w%failing_call, int32), 2, "close is the failing call")
        call assert_no_file(filename, "no root leaves no file")
    end subroutine test_no_root_at_close

    !> The first error stays; later calls, even ones that would fail differently, change nothing.
    subroutine test_errors_are_sticky()
        character(len=*), parameter :: filename = "json_seq_sticky.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call open_fresh(w, filename)
        call json_begin_array(w)
        call json_member(w, "k", 1_int32)
        call json_element(w, ieee_value(1.0_real64, ieee_quiet_nan))
        call json_element(w, "ok"//char(255))
        call json_end_object(w)
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "the first error is returned")
        call assert_equal_int(int(w%failing_call, int32), 3, "the first failing call is kept")
        call assert_no_file(filename, "sticky failure leaves no file")
    end subroutine test_errors_are_sticky

    subroutine test_writer_not_open()
        type(json_writer) :: unopened_for_value, unopened_for_close
        integer(int32) :: ierr

        call json_element(unopened_for_value, 1_int32)
        call assert_err(unopened_for_value%ierr, ERR_JSON_SEQUENCE, "value on a writer never opened")

        call json_close(unopened_for_close, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "close on a writer never opened")
    end subroutine test_writer_not_open

    !> Opening an open writer fails it and removes the file it was writing.
    subroutine test_open_twice()
        character(len=*), parameter :: first_file = "json_open_twice_1.test.json"
        character(len=*), parameter :: second_file = "json_open_twice_2.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call remove_file(first_file)
        call remove_file(second_file)
        call json_open(w, first_file, ierr)
        call json_begin_array(w)
        call json_open(w, second_file, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "second open")
        call assert_equal_int(int(w%failing_call, int32), 3, "the second open is the failing call")
        call assert_no_file(first_file, "the first file is deleted")
        call assert_no_file(second_file, "the second file is never created")
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "close returns the first error")
    end subroutine test_open_twice

    !> A closed writer, failed or not, may be opened again and starts from scratch.
    subroutine test_writer_reuse()
        character(len=*), parameter :: filename = "json_reuse.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call open_fresh(w, filename)
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, "first use fails")

        call json_open(w, filename, ierr)
        call assert_err(ierr, ERR_OK, "reopen after a failure")
        call assert_equal_int(int(w%failing_call, int32), 0, "the failure is forgotten")
        call json_element(w, 1_int32)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "second use succeeds")
        call assert_file_equals(filename, "1"//LF, "second document")
        call remove_file(filename)

        call json_open(w, filename, ierr)
        call json_element(w, 2_int32)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "third use succeeds")
        call assert_file_equals(filename, "2"//LF, "third document")

        call json_element(w, 3_int32)
        call assert_err(w%ierr, ERR_JSON_SEQUENCE, "a value after close")
        call assert_file_equals(filename, "2"//LF, "a call after close leaves the finished file alone")
        call remove_file(filename)
    end subroutine test_writer_reuse

    ! ============================================================================================
    ! The file
    ! ============================================================================================

    subroutine test_existing_file_refused()
        character(len=*), parameter :: filename = "json_existing.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, unit

        call remove_file(filename)
        open (newunit=unit, file=filename, access="stream", form="unformatted", status="new", action="write")
        write (unit) "keep me"
        close (unit)

        call json_open(w, filename, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "an existing file is refused")
        call assert_equal_int(int(w%failing_call, int32), 1, "open is the failing call")
        call json_element(w, 1_int32)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "close returns the open error")
        call assert_file_equals(filename, "keep me", "the existing file is untouched")
        call remove_file(filename)
    end subroutine test_existing_file_refused

    subroutine test_unwritable_path()
        character(len=*), parameter :: filename = "json_no_such_directory.test/out.json"
        type(json_writer) :: w
        integer(int32) :: ierr

        call json_open(w, filename, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "a path in a missing directory")
        call json_close(w, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "close returns the open error")
        call assert_no_file(filename, "nothing is created")
    end subroutine test_unwritable_path

    !> Nesting deeper than the initial stack grows it.
    subroutine test_deep_nesting()
        character(len=*), parameter :: filename = "json_deep.test.json"
        integer(int32), parameter :: depth = 100
        type(json_writer) :: w
        integer(int32) :: ierr, level

        call remove_file(filename)
        call json_open(w, filename, ierr)
        do level = 1, depth
            if (mod(level, 2) == 1) then
                call json_begin_array(w)
            else
                call json_begin_object(w)
                call json_begin_array(w, "k")
                call json_end_array(w)
                call json_begin_array(w, "a")
            end if
        end do
        call json_element(w, 1_int32)
        do level = depth, 1, -1
            call json_end_array(w)
            if (mod(level, 2) == 0) call json_end_object(w)
        end do
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "deep nesting")
        call assert_file_equals(filename, repeat('[{"k":[],"a":[', depth/2)//"1"//repeat("]}]", depth/2)//LF, &
                                "deep nesting")
        call remove_file(filename)
    end subroutine test_deep_nesting

    !> A string larger than the whole buffer bypasses it.
    subroutine test_string_longer_than_buffer()
        character(len=*), parameter :: filename = "json_long_string.test.json"
        integer(int32), parameter :: text_length = 300000
        type(json_writer) :: w
        integer(int32) :: ierr
        character(len=:), allocatable :: text

        allocate (character(len=text_length) :: text)
        text(:) = repeat("abcdefghij", text_length/10)
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, "head")
        call json_element(w, text)
        call json_element(w, "tail")
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "long string")
        call assert_file_equals(filename, '["head","'//text//'","tail"]'//LF, "long string")
        call remove_file(filename)
    end subroutine test_string_longer_than_buffer

    !> A document of several buffers' worth is written in full.
    subroutine test_output_larger_than_buffer()
        character(len=*), parameter :: filename = "json_large.test.json"
        integer(int32), parameter :: n_values = 40000
        type(json_writer) :: w
        integer(int32) :: ierr, i_value
        real(real64), allocatable :: values(:), parsed(:)

        allocate (values(n_values), parsed(n_values))
        do i_value = 1, n_values
            values(i_value) = real(i_value, real64)/7.0_real64
        end do
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_element(w, values)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "large document")
        call read_real_array_file(filename, parsed, ierr)
        call assert_equal_int(ierr, 0, "large document parses")
        call assert_equal_int(count_bit_mismatches(parsed, values), 0, "large document reads back bit for bit")
        call remove_file(filename)
    end subroutine test_output_larger_than_buffer

    !> c_bool with the byte value C writes for true (1), and default logical arrays.
    subroutine test_boolean_kinds()
        character(len=*), parameter :: filename = "json_booleans.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr
        logical(c_bool) :: c_true(2)

        c_true(1) = transfer(1_int8, c_true(1))
        c_true(2) = transfer(0_int8, c_true(2))
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, c_true(1))
        call json_element(w, c_true)
        call json_element(w, [.true., .false.])
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "booleans")
        call assert_file_equals(filename, "[true,[true,false],[true,false]]"//LF, "booleans of both kinds")
        call remove_file(filename)
    end subroutine test_boolean_kinds

    !> Values of a type the writer does not support fail it, as scalars and as arrays, keyed or not.
    subroutine test_unsupported_types()
        character(len=*), parameter :: filename = "json_unsupported.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, i_case
        real(real32) :: singles(2)
        integer(int16) :: shorts(3)

        singles = 1.0_real32
        shorts = 2_int16
        do i_case = 1, 6
            call remove_file(filename)
            call json_open(w, filename, ierr)
            call json_begin_object(w)
            call json_member(w, "ok", 1_int32)
            select case (i_case)
            case (1)
                call json_member(w, "x", 1.0_real32)
            case (2)
                call json_member(w, "x", 1_int8)
            case (3)
                call json_member(w, "x", (1.0_real64, 2.0_real64))
            case (4)
                call json_member(w, "x", singles)
            case (5)
                call json_member(w, "x", shorts)
            case (6)
                call json_begin_array(w, "a")
                call json_element(w, shorts(2:3))
            end select
            call json_close(w, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "unsupported type, case "//i_case)
            call assert_equal_int(int(w%failing_call, int32), merge(5, 4, i_case == 6), &
                                  "the call with the unsupported value fails, case "//i_case)
            call assert_no_file(filename, "unsupported type leaves no file, case "//i_case)
        end do
    end subroutine test_unsupported_types

    !> Zeros are recognised from the bits: a negative zero keeps its sign, and the smallest
    !| subnormals, which compare equal to zero where denormals are flushed, are written in full.
    subroutine test_zero_decided_from_bits()
        character(len=*), parameter :: filename = "json_zero_bits.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr
        real(real64) :: smallest_subnormal

        smallest_subnormal = transfer(1_int64, smallest_subnormal)
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, "negative_zero", -0.0_real64)
        call json_member(w, "subnormal", smallest_subnormal)
        call json_member(w, "negative_subnormal", -smallest_subnormal)
        call json_member(w, "array", [smallest_subnormal, -0.0_real64, 0.0_real64])
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "zeros and subnormals")
        call assert_file_equals(filename, '{"negative_zero":-0.0,"subnormal":4.9406564584124654E-324,' &
                                //'"negative_subnormal":-4.9406564584124654E-324,' &
                                //'"array":[4.9406564584124654E-324,-0.0,0.0]}'//LF, "zeros and subnormals")
        call remove_file(filename)
    end subroutine test_zero_decided_from_bits

    !> An empty array is [] whatever its element type, a type the writer cannot write included.
    subroutine test_empty_arrays()
        character(len=*), parameter :: filename = "json_empty_arrays.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr
        integer(int32) :: values(5)

        values = 1
        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_object(w)
        call json_member(w, "int32", [integer(int32) ::])
        call json_member(w, "int64", [integer(int64) ::])
        call json_member(w, "real64", [real(real64) ::])
        call json_member(w, "c_bool", [logical(c_bool) ::])
        call json_member(w, "logical", [logical ::])
        call json_member(w, "character", [character(len=5) ::])
        call json_member(w, "section", values(4:3))
        call json_member(w, "real32", [real(real32) ::])
        call json_begin_array(w, "elements")
        call json_element(w, [integer ::])
        call json_element(w, values(5:1))
        call json_end_array(w)
        call json_end_object(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "empty arrays")
        call assert_file_equals(filename, '{"int32":[],"int64":[],"real64":[],"c_bool":[],"logical":[],' &
                                //'"character":[],"section":[],"real32":[],"elements":[[],[]]}'//LF, "empty arrays")
        call remove_file(filename)
    end subroutine test_empty_arrays

    !> Sections that are not one contiguous run: stride 2, a negative stride, a non-unit lower
    !| bound, a matrix row and a matrix column, for integers, reals, booleans and strings. The
    !| expected text is built by walking the indices explicitly.
    subroutine test_array_sections()
        character(len=*), parameter :: filename = "json_sections.test.json"
        type(json_writer) :: w
        integer(int32) :: ierr, i, j
        integer(int32) :: numbers(-7:2), number_matrix(3, 4)
        real(real64) :: reals(-7:2), real_matrix(3, 4)
        logical(c_bool) :: flags(-7:2), flag_matrix(3, 4)
        character(len=2) :: letters(-7:2), letter_matrix(3, 4)
        character(len=:), allocatable :: expected

        do i = -7, 2
            numbers(i) = 10*i
            reals(i) = real(i, real64)
            flags(i) = mod(i + 9, 3) == 1
            letters(i) = achar(iachar("a") + i + 7)
        end do
        do j = 1, 4
            do i = 1, 3
                number_matrix(i, j) = 10*i + j
                real_matrix(i, j) = real(i + 3*(j - 1) - 6, real64)
                flag_matrix(i, j) = mod(i + j, 2) == 0
                letter_matrix(i, j) = achar(iachar("A") + i - 1)//achar(iachar("0") + j)
            end do
        end do

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call json_begin_array(w)
        call json_element(w, numbers(-7:2:2))
        call json_element(w, numbers(2:-7:-3))
        call json_element(w, numbers(-5:-3))
        call json_element(w, numbers)
        call json_element(w, number_matrix(2, :))
        call json_element(w, number_matrix(:, 3))
        call json_element(w, reals(-7:2:2))
        call json_element(w, reals(2:-7:-3))
        call json_element(w, real_matrix(2, :))
        call json_element(w, real_matrix(:, 3))
        call json_element(w, flags(-7:2:2))
        call json_element(w, flags(2:-7:-3))
        call json_element(w, flag_matrix(2, :))
        call json_element(w, flag_matrix(:, 3))
        call json_element(w, letters(-7:2:2))
        call json_element(w, letters(2:-7:-3))
        call json_element(w, letter_matrix(2, :))
        call json_element(w, letter_matrix(:, 3))
        call json_end_array(w)
        call json_close(w, ierr)
        call assert_err(ierr, ERR_OK, "sections")

        expected = "["//int_list(numbers, -7, -7, 2, 2)//","//int_list(numbers, -7, 2, -7, -3)//"," &
                   //int_list(numbers, -7, -5, -3, 1)//","//int_list(numbers, -7, -7, 2, 1)//"," &
                   //int_list(number_matrix(2, :), 1, 1, 4, 1)//","//int_list(number_matrix(:, 3), 1, 1, 3, 1)
        expected = expected//","//real_list(reals, -7, -7, 2, 2)//","//real_list(reals, -7, 2, -7, -3)//"," &
                   //real_list(real_matrix(2, :), 1, 1, 4, 1)//","//real_list(real_matrix(:, 3), 1, 1, 3, 1)
        expected = expected//","//flag_list(flags, -7, -7, 2, 2)//","//flag_list(flags, -7, 2, -7, -3)//"," &
                   //flag_list(flag_matrix(2, :), 1, 1, 4, 1)//","//flag_list(flag_matrix(:, 3), 1, 1, 3, 1)
        expected = expected//","//letter_list(letters, -7, -7, 2, 2)//","//letter_list(letters, -7, 2, -7, -3)//"," &
                   //letter_list(letter_matrix(2, :), 1, 1, 4, 1)//","//letter_list(letter_matrix(:, 3), 1, 1, 3, 1)//"]"//LF
        call assert_file_equals(filename, expected, "sections")
        ! a spot check of the builders themselves, so the comparison cannot pass by agreement alone
        call assert_string_equal(int_list(numbers, -7, 2, -7, -3), "[20,-10,-40,-70]", "negative stride, by hand")
        call assert_string_equal(letter_list(letter_matrix(2, :), 1, 1, 4, 1), '["B1","B2","B3","B4"]', "row, by hand")
        call remove_file(filename)
    end subroutine test_array_sections

    !> `[values(first), values(first+step), ...]` for integers, walking the indices one by one;
    !| `lower` is the index the first element of `values` has in the caller.
    function int_list(values, lower, first, last, step) result(text)
        integer(int32), intent(in) :: values(:)
        integer, intent(in) :: lower, first, last, step
        character(len=:), allocatable :: text
        character(len=16) :: digits
        integer :: i

        text = "["
        do i = first, last, step
            write (digits, "(I0)") values(i - lower + 1)
            if (i /= first) text = text//","
            text = text//trim(digits)
        end do
        text = text//"]"
    end function int_list

    !> The same for reals holding small whole numbers, whose 17-digit text is written out here.
    function real_list(values, lower, first, last, step) result(text)
        real(real64), intent(in) :: values(:)
        integer, intent(in) :: lower, first, last, step
        character(len=:), allocatable :: text
        integer :: i, whole

        text = "["
        do i = first, last, step
            whole = nint(values(i - lower + 1))
            if (i /= first) text = text//","
            if (whole == 0) then
                text = text//"0.0"
            else
                ! a digit by substring: nvfortran 26.9 pads achar() of a variable in a concatenation
                if (whole < 0) text = text//"-"
                text = text//"0123456789"(abs(whole) + 1:abs(whole) + 1)//".0000000000000000E+000"
            end if
        end do
        text = text//"]"
    end function real_list

    function flag_list(values, lower, first, last, step) result(text)
        logical(c_bool), intent(in) :: values(:)
        integer, intent(in) :: lower, first, last, step
        character(len=:), allocatable :: text
        integer :: i

        text = "["
        do i = first, last, step
            if (i /= first) text = text//","
            if (values(i - lower + 1)) then
                text = text//"true"
            else
                text = text//"false"
            end if
        end do
        text = text//"]"
    end function flag_list

    function letter_list(values, lower, first, last, step) result(text)
        character(len=*), intent(in) :: values(:)
        integer, intent(in) :: lower, first, last, step
        character(len=:), allocatable :: text
        integer :: i

        text = "["
        do i = first, last, step
            if (i /= first) text = text//","
            text = text//'"'//trim(values(i - lower + 1))//'"'
        end do
        text = text//"]"
    end function letter_list


    ! ============================================================================================
    ! Random call sequences against an independent model of the rules
    ! ============================================================================================

    !> Seeded random call sequences, about 3% of them with a wrong key choice, some with NaN or
    !| invalid UTF-8. A small model of the rules predicts the text, or the first error and its
    !| call number; the writer must agree on every sequence.
    subroutine test_fuzz_call_sequences()
        character(len=*), parameter :: filename = "json_fuzz.test.json"
        integer(int32), parameter :: n_sequences = 400, max_calls = 40, max_depth = 64
        integer(int32), parameter :: OBJECT = 1, ARRAY = 2
        type(json_writer) :: w
        integer(int64) :: random_state
        integer(int32) :: i_sequence, i_call, n_calls, operation, ierr, i_pick, expected_ierr, expected_call
        integer(int32) :: depth, call_number, n_mismatches
        integer(int32) :: kinds(max_depth)
        logical :: with_key, is_first, root_done, is_value, closes_ok
        character(len=:), allocatable :: expected, piece, got
        logical :: found
        character(len=8) :: key_pool(3), key_text(3)
        character(len=12) :: string_pool(4)
        character(len=16) :: string_text(4)

        ! the backslash is spelled achar(92), which no compiler reads as an escape
        key_pool = [character(len=8) :: "k", 'q"', ""]
        key_text = [character(len=8) :: '"k"', '"q'//BS//'""', '""']
        string_pool = [character(len=12) :: "plain", "tab"//achar(9), "back"//BS//"slash", ""]
        string_text = [character(len=16) :: '"plain"', '"tab'//BS//'t"', '"back'//BS//BS//'slash"', '""']

        random_state = 20260923_int64
        n_mismatches = 0
        do i_sequence = 1, n_sequences
            call remove_file(filename)
            call json_open(w, filename, ierr)
            expected = ""
            expected_ierr = ERR_OK
            expected_call = 0
            depth = 0
            is_first = .true.
            root_done = .false.
            call_number = 1
            n_calls = next_random(random_state, max_calls + 1) - 1
            do i_call = 1, n_calls
                call_number = call_number + 1
                operation = next_random(random_state, 10)
                ! the right key choice, flipped for about 3% of the calls
                with_key = depth > 0
                if (with_key) with_key = kinds(depth) == OBJECT
                if (next_random(random_state, 100) <= 3) with_key = .not. with_key
                i_pick = next_random(random_state, 3)
                if (depth >= max_depth .and. operation <= 2) operation = 5
                is_value = operation /= 3 .and. operation /= 4

                piece = ""
                select case (operation)
                case (1)
                    if (with_key) then
                        call json_begin_object(w, trim(key_pool(i_pick)))
                    else
                        call json_begin_object(w)
                    end if
                    piece = "{"
                case (2)
                    if (with_key) then
                        call json_begin_array(w, trim(key_pool(i_pick)))
                    else
                        call json_begin_array(w)
                    end if
                    piece = "["
                case (3)
                    call json_end_object(w)
                case (4)
                    call json_end_array(w)
                case (5)
                    if (with_key) then
                        call json_member(w, trim(key_pool(i_pick)), -12345_int32)
                    else
                        call json_element(w, -12345_int32)
                    end if
                    piece = "-12345"
                case (6)
                    if (with_key) then
                        call json_member(w, trim(key_pool(i_pick)), string_pool(i_pick + 1))
                    else
                        call json_element(w, string_pool(i_pick + 1))
                    end if
                    piece = trim(string_text(i_pick + 1))
                case (7)
                    if (with_key) then
                        call json_member(w, trim(key_pool(i_pick)), .true._c_bool)
                    else
                        call json_element(w, .true._c_bool)
                    end if
                    piece = "true"
                case (8)
                    if (with_key) then
                        call json_null(w, trim(key_pool(i_pick)))
                    else
                        call json_null(w)
                    end if
                    piece = "null"
                case (9)
                    if (with_key) then
                        call json_member(w, trim(key_pool(i_pick)), [0.5_real64, -0.0_real64])
                    else
                        call json_element(w, [0.5_real64, -0.0_real64])
                    end if
                    piece = "[5.0000000000000000E-001,-0.0]"
                case (10)
                    ! rarely a value the writer must refuse for its content
                    if (next_random(random_state, 10) == 1) then
                        if (with_key) then
                            call json_member(w, trim(key_pool(i_pick)), ieee_value(1.0_real64, ieee_quiet_nan))
                        else
                            call json_element(w, ieee_value(1.0_real64, ieee_quiet_nan))
                        end if
                        piece = "NaN"
                    else if (next_random(random_state, 10) == 1) then
                        if (with_key) then
                            call json_member(w, trim(key_pool(i_pick)), "x"//char(192)//char(128))
                        else
                            call json_element(w, "x"//char(192)//char(128))
                        end if
                        piece = "UTF8"
                    else
                        if (with_key) then
                            call json_member(w, trim(key_pool(i_pick)), 7_int64)
                        else
                            call json_element(w, 7_int64)
                        end if
                        piece = "7"
                    end if
                end select

                ! the model
                if (expected_ierr /= ERR_OK) cycle
                if (is_value) then
                    closes_ok = .true.
                    if (depth == 0) then
                        if (root_done .or. with_key) closes_ok = .false.
                    else
                        if ((kinds(depth) == OBJECT) .neqv. with_key) closes_ok = .false.
                    end if
                    if (.not. closes_ok) then
                        expected_ierr = ERR_JSON_SEQUENCE
                    else if (piece == "NaN") then
                        expected_ierr = ERR_NAN_INF
                    else if (piece == "UTF8") then
                        expected_ierr = ERR_INVALID_UTF8
                    end if
                    if (expected_ierr /= ERR_OK) then
                        expected_call = call_number
                        cycle
                    end if
                    if (.not. is_first) expected = expected//","
                    if (with_key) expected = expected//trim(key_text(i_pick))//":"
                    expected = expected//piece
                    is_first = .false.
                    if (operation <= 2) then
                        depth = depth + 1
                        kinds(depth) = merge(OBJECT, ARRAY, operation == 1)
                        is_first = .true.
                    else if (depth == 0) then
                        root_done = .true.
                    end if
                else
                    closes_ok = depth > 0
                    if (closes_ok) closes_ok = kinds(depth) == merge(OBJECT, ARRAY, operation == 3)
                    if (.not. closes_ok) then
                        expected_ierr = ERR_JSON_SEQUENCE
                        expected_call = call_number
                        cycle
                    end if
                    expected = expected//merge("}", "]", operation == 3)
                    depth = depth - 1
                    is_first = .false.
                    if (depth == 0) root_done = .true.
                end if
            end do
            call_number = call_number + 1
            if (expected_ierr == ERR_OK .and. (depth > 0 .or. .not. root_done)) then
                expected_ierr = ERR_JSON_SEQUENCE
                expected_call = call_number
            end if

            call json_close(w, ierr)
            if (ierr /= expected_ierr .or. w%failing_call /= expected_call) then
                n_mismatches = n_mismatches + 1
                call assert_equal_int(ierr, expected_ierr, "fuzz sequence "//i_sequence//": error")
                call assert_equal_int(int(w%failing_call, int32), expected_call, "fuzz sequence "//i_sequence//": call")
            end if
            if (expected_ierr == ERR_OK) then
                call read_whole_file(filename, got, found)
                if (.not. found .or. got /= expected//LF .or. len(got) /= len(expected) + 1) then
                    n_mismatches = n_mismatches + 1
                    call assert_file_equals(filename, expected//LF, "fuzz sequence "//i_sequence//": text")
                end if
            else
                if (file_exists(filename)) n_mismatches = n_mismatches + 1
                call assert_no_file(filename, "fuzz sequence "//i_sequence//": no file")
            end if
            if (n_mismatches > 5) exit
        end do
        call remove_file(filename)
    end subroutine test_fuzz_call_sequences

    ! ============================================================================================
    ! Helpers
    ! ============================================================================================

    !> A random integer from 1 to `n_choices`, by the Park-Miller generator (no overflow in int64).
    function next_random(state, n_choices) result(choice)
        integer(int64), intent(inout) :: state
        integer(int32), intent(in) :: n_choices
        integer(int32) :: choice

        state = mod(state*48271_int64, 2147483647_int64)
        choice = int(mod(state, int(n_choices, int64)), int32) + 1
    end function next_random

    !> How many of two equally long arrays differ in their bits.
    integer(int32) function count_bit_mismatches(actual, expected)
        real(real64), intent(in) :: actual(:), expected(:)
        integer :: i_value

        count_bit_mismatches = 0
        do i_value = 1, size(expected)
            if (transfer(actual(i_value), 0_int64) /= transfer(expected(i_value), 0_int64)) &
                count_bit_mismatches = count_bit_mismatches + 1
        end do
    end function count_bit_mismatches

    !> The string made of the given byte values.
    function bytes(values) result(text)
        integer, intent(in) :: values(:)
        character(len=size(values)) :: text
        integer :: i_byte

        do i_byte = 1, size(values)
            text(i_byte:i_byte) = achar(values(i_byte))
        end do
    end function bytes

    integer(int32) function n_invalid_utf8_cases()
        n_invalid_utf8_cases = 21
    end function n_invalid_utf8_cases

    !> One malformed UTF-8 sequence per case: stray continuations, overlongs, surrogates, beyond
    !| U+10FFFF, truncations, and a continuation that is not one.
    function invalid_utf8_case(i_case) result(text)
        integer(int32), intent(in) :: i_case
        character(len=8) :: text

        text = ""
        select case (i_case)
        case (1); text = bytes([128])                 ! lone continuation
        case (2); text = bytes([191])                 ! lone continuation
        case (3); text = bytes([192, 128])            ! overlong NUL
        case (4); text = bytes([193, 191])            ! overlong U+007F
        case (5); text = bytes([224, 128, 128])       ! overlong 3-byte
        case (6); text = bytes([224, 159, 191])       ! overlong U+07FF
        case (7); text = bytes([240, 128, 128, 128])  ! overlong 4-byte
        case (8); text = bytes([240, 143, 191, 191])  ! overlong U+FFFF
        case (9); text = bytes([237, 160, 128])       ! surrogate U+D800
        case (10); text = bytes([237, 191, 191])      ! surrogate U+DFFF
        case (11); text = bytes([244, 144, 128, 128]) ! U+110000
        case (12); text = bytes([245, 128, 128, 128]) ! lead byte beyond U+10FFFF
        case (13); text = bytes([255])                ! never valid
        case (14); text = bytes([195])                ! truncated 2-byte
        case (15); text = bytes([226, 130])           ! truncated 3-byte
        case (16); text = bytes([240, 159, 152])      ! truncated 4-byte
        case (17); text = bytes([195, 65])            ! continuation that is ASCII
        case (18); text = bytes([226, 130, 32, 172])  ! blank inside a sequence
        case (19); text = bytes([226, 130, 65])       ! ASCII as the third byte of three
        case (20); text = bytes([240, 159, 152, 65])  ! ASCII as the fourth byte of four
        case (21); text = bytes([240, 159, 65, 128])  ! ASCII as the third byte of four
        end select
    end function invalid_utf8_case

    !> Remove a leftover file and open the writer on it.
    subroutine open_fresh(w, filename)
        type(json_writer), intent(inout) :: w
        character(len=*), intent(in) :: filename
        integer(int32) :: ierr

        call remove_file(filename)
        call json_open(w, filename, ierr)
        call assert_err(ierr, ERR_OK, "open "//filename)
    end subroutine open_fresh

    !> After a sequence violation: the error, the call that caused it, no file, and close agrees.
    subroutine assert_sequence_failure(w, filename, failing_call, msg)
        type(json_writer), intent(inout) :: w
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: failing_call
        character(len=*), intent(in) :: msg
        integer(int32) :: ierr

        call assert_err(w%ierr, ERR_JSON_SEQUENCE, msg)
        call assert_equal_int(int(w%failing_call, int32), failing_call, msg//": failing call")
        call assert_no_file(filename, msg//": file deleted at once")
        call json_close(w, ierr)
        call assert_err(ierr, ERR_JSON_SEQUENCE, msg//": close returns it")
        call assert_equal_int(int(w%failing_call, int32), failing_call, msg//": close keeps the call")
        call assert_no_file(filename, msg//": no file after close")
    end subroutine assert_sequence_failure

    !> Read a file holding one JSON array of numbers back into reals.
    subroutine read_real_array_file(filename, values, ierr)
        character(len=*), intent(in) :: filename
        real(real64), intent(out) :: values(:)
        integer(int32), intent(out) :: ierr
        character(len=:), allocatable :: content
        logical :: found
        integer :: array_end

        values = 0.0_real64
        call read_whole_file(filename, content, found)
        ierr = 1
        if (.not. found) return
        array_end = index(content, "]", back=.true.)
        if (content(1:1) /= "[" .or. array_end < 2) return
        ! list-directed input takes the commas as separators
        read (content(2:array_end - 1), *, iostat=ierr) values
    end subroutine read_real_array_file

    !> The whole content of a file, byte for byte.
    subroutine read_whole_file(filename, content, found)
        character(len=*), intent(in) :: filename
            !! Path of the file
        character(len=:), allocatable, intent(out) :: content
            !! Its bytes; empty when it does not exist
        logical, intent(out) :: found
            !! Whether it exists
        integer(int64) :: file_size
        integer :: unit, iostat

        inquire (file=filename, exist=found)
        if (.not. found) then
            allocate (character(len=0) :: content)
            return
        end if
        inquire (file=filename, size=file_size)
        allocate (character(len=file_size) :: content)
        open (newunit=unit, file=filename, access="stream", form="unformatted", status="old", action="read", &
              iostat=iostat)
        if (iostat /= 0) return
        if (file_size > 0) read (unit, iostat=iostat) content
        close (unit)
    end subroutine read_whole_file

    logical function file_exists(filename)
        character(len=*), intent(in) :: filename
        inquire (file=filename, exist=file_exists)
    end function file_exists

    subroutine remove_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, iostat

        if (.not. file_exists(filename)) return
        open (newunit=unit, file=filename, status="old", iostat=iostat)
        if (iostat == 0) close (unit, status="delete")
    end subroutine remove_file

    !> Assert that a file holds exactly `expected`, showing the first difference otherwise.
    subroutine assert_file_equals(filename, expected, msg)
        character(len=*), intent(in) :: filename, expected, msg
        character(len=:), allocatable :: content
        logical :: found
        integer :: first_difference, shown_end

        call read_whole_file(filename, content, found)
        if (.not. found) then
            call assert_true(.false., msg//": file '"//filename//"' does not exist")
            return
        end if
        if (len(content) == len(expected)) then
            if (content == expected) return
        end if
        first_difference = 1
        do while (first_difference <= min(len(content), len(expected)))
            if (content(first_difference:first_difference) /= expected(first_difference:first_difference)) exit
            first_difference = first_difference + 1
        end do
        shown_end = min(len(content), first_difference + 60)
        call assert_true(.false., msg//": content differs from byte "//first_difference//" (length " &
                         //len(content)//", expected "//len(expected)//"): got '" &
                         //content(max(1, first_difference - 20):shown_end)//"', expected '" &
                         //expected(max(1, first_difference - 20):min(len(expected), first_difference + 60))//"'")
    end subroutine assert_file_equals

    subroutine assert_no_file(filename, msg)
        character(len=*), intent(in) :: filename, msg
        call assert_false(file_exists(filename), msg//": '"//filename//"' must not exist")
    end subroutine assert_no_file

end module mod_test_json_serialize
