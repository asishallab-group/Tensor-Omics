!> What the children of the `tox_data_read_write` suite share: the published file layout, written
!| out by hand, and the checks every `save_*` and `load_*` procedure gets, once per kind of array.
!|
!| A file is an unformatted stream of int32 words and elements: the magic number 1178677808 (hex
!| 46413230, 'FA20' read as one word), the type code (-1 integer, -2 real, a string length for
!| strings), the number of dimensions, the extents, then the elements in column-major order
!| (f42_serde_arrays_utils documents it). The words are in the machine's byte order -- on a
!| little-endian machine the magic number's bytes read "02AF" -- and so is every file this suite
!| writes or expects: the layout as pinned here is the native one, not a portable one.
!|
!| `write_*_file` write that layout from the I/O list alone, independently of the serde code.
!| Each `save_*` case compares a written file byte for byte with the file `write_*_file` makes of
!| the same array, and each `load_*` case reads a file `write_*_file` made, so neither direction
!| is tested only through the other.
!|
!| The procedures of one kind take the same arguments, so each check takes the procedure under
!| test as an argument: `tox_data_read_write`'s three real matrices share one set of checks, its
!| two string arrays another. Each check is one case of a child, in a scratch directory of its
!| own, named after the case.
!|
!| A missing file and a path that cannot be written are refused with ERR_FILE_OPEN, and the error
!| unit stays connected (the serde layer used to close the unit number it had not opened).
module mod_test_tox_data_read_write_fixtures
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int8
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_errors
    use scratch_directory
    implicit none
    private
    public :: save_strings_interface, load_strings_interface, save_integers_interface, load_integers_interface
    public :: save_reals_interface, load_reals_interface
    public :: check_saves_strings, check_save_replaces_strings, check_loads_strings, check_load_refuses_strings
    public :: check_round_trip_strings
    public :: check_saves_integers, check_save_replaces_integers, check_loads_integers, check_load_refuses_integers
    public :: check_round_trip_integers
    public :: check_saves_reals, check_save_replaces_reals, check_loads_reals, check_load_refuses_reals
    public :: check_round_trip_reals, check_saves_non_finite_reals, check_loads_non_finite_reals
    public :: check_save_refuses_an_unwritable_path_strings, check_load_refuses_a_missing_file_strings
    public :: check_save_refuses_an_unwritable_path_integers, check_load_refuses_a_missing_file_integers
    public :: check_save_refuses_an_unwritable_path_reals, check_load_refuses_a_missing_file_reals

    !> The magic number every array file starts with, 'FA20' read as an int32.
    integer(int32), parameter :: MAGIC = 1178677808
    integer(int32), parameter :: INTEGER_CODE = -1, REAL_CODE = -2

    !> The broken files every `load_*` procedure refuses, how each is broken, and its code.
    integer(int32), parameter :: N_BAD_FILES = 12
    character(len=32), parameter :: BAD_FILE_NAMES(N_BAD_FILES) = [character(len=32) :: &
        "empty", "part of the magic number", "a wrong magic number", "only the magic number", &
        "no number of dimensions", "16 dimensions", "-1 dimensions", "an extent missing", &
        "another type", "another rank", "another extent", "an element short"]
    integer(int32), parameter :: BAD_FILE_CODES(N_BAD_FILES) = [ &
        ERR_READ_MAGIC, ERR_READ_MAGIC, ERR_INVALID_FORMAT, ERR_READ_TYPE, &
        ERR_READ_NDIMS, ERR_INVALID_FORMAT, ERR_INVALID_FORMAT, ERR_READ_DIMS, &
        ERR_TYPE_MISMATCH, ERR_DIM_MISMATCH, ERR_DIM_MISMATCH, ERR_READ_DATA]

    abstract interface
        subroutine save_strings_interface(strings, filename, ierr)
            import :: int32
            character(len=*), intent(in), contiguous :: strings(:)
            character(len=*), intent(in) :: filename
            integer(int32), intent(out) :: ierr
        end subroutine save_strings_interface

        subroutine load_strings_interface(strings, filename, ierr)
            import :: int32
            character(len=*), intent(in) :: filename
            character(len=*), intent(out), contiguous :: strings(:)
            integer(int32), intent(out) :: ierr
        end subroutine load_strings_interface

        subroutine save_integers_interface(integers, filename, ierr)
            import :: int32
            character(len=*), intent(in) :: filename
            integer(int32), intent(in), contiguous :: integers(:)
            integer(int32), intent(out) :: ierr
        end subroutine save_integers_interface

        subroutine load_integers_interface(integers, filename, ierr)
            import :: int32
            character(len=*), intent(in) :: filename
            integer(int32), intent(out), contiguous :: integers(:)
            integer(int32), intent(out) :: ierr
        end subroutine load_integers_interface

        subroutine save_reals_interface(reals, filename, ierr)
            import :: int32, real64
            character(len=*), intent(in) :: filename
            real(real64), intent(in), contiguous :: reals(:, :)
            integer(int32), intent(out) :: ierr
        end subroutine save_reals_interface

        subroutine load_reals_interface(reals, filename, ierr)
            import :: int32, real64
            character(len=*), intent(in) :: filename
            real(real64), intent(out), contiguous :: reals(:, :)
            integer(int32), intent(out) :: ierr
        end subroutine load_reals_interface
    end interface

contains

    ! ---- the layout, written from the I/O list ----

    !> Opens `file_path` for a stream to replace it, failing an assertion if it cannot.
    subroutine open_for_writing(file_path, unit, status)
        character(*), intent(in) :: file_path
        integer(int32), intent(out) :: unit, status

        open (newunit=unit, file=file_path, access='stream', form='unformatted', status='replace', action='write', &
              iostat=status)
        call assert_equal_int(status, 0_int32, "open_for_writing: "//file_path)
    end subroutine open_for_writing

    !> The layout of a 1-d array of strings: the type code is their length.
    subroutine write_strings_file(file_path, strings)
        character(*), intent(in) :: file_path
        character(len=*), intent(in) :: strings(:)
        integer(int32) :: unit, status

        call open_for_writing(file_path, unit, status)
        if (status /= 0) return
        write (unit) MAGIC, len(strings, kind=int32), 1_int32, size(strings, kind=int32)
        write (unit) strings
        close (unit)
    end subroutine write_strings_file

    !> The layout of a 1-d int32 array.
    subroutine write_integers_file(file_path, integers)
        character(*), intent(in) :: file_path
        integer(int32), intent(in) :: integers(:)
        integer(int32) :: unit, status

        call open_for_writing(file_path, unit, status)
        if (status /= 0) return
        write (unit) MAGIC, INTEGER_CODE, 1_int32, size(integers, kind=int32)
        write (unit) integers
        close (unit)
    end subroutine write_integers_file

    !> The layout of a 2-d real64 array.
    subroutine write_reals_file(file_path, reals)
        character(*), intent(in) :: file_path
        real(real64), intent(in), contiguous :: reals(:, :)
        integer(int32) :: unit, status

        call open_for_writing(file_path, unit, status)
        if (status /= 0) return
        write (unit) MAGIC, REAL_CODE, 2_int32, size(reals, 1, kind=int32), size(reals, 2, kind=int32)
        write (unit) reals
        close (unit)
    end subroutine write_reals_file

    !> All bytes of `file_path`; none for a file that cannot be read.
    subroutine read_file_bytes(file_path, bytes)
        character(*), intent(in) :: file_path
        integer(int8), allocatable, intent(out) :: bytes(:)
        integer(int32) :: unit, status, file_size

        open (newunit=unit, file=file_path, access='stream', form='unformatted', status='old', action='read', &
              iostat=status)
        if (status /= 0) then
            allocate (bytes(0))
            return
        end if
        inquire (unit=unit, size=file_size)
        allocate (bytes(file_size))
        if (file_size > 0) read (unit, iostat=status) bytes
        close (unit)
    end subroutine read_file_bytes

    !> Fails unless `actual_path` holds exactly the bytes of `expected_path`.
    subroutine check_same_file(actual_path, expected_path, msg)
        character(*), intent(in) :: actual_path, expected_path, msg
        integer(int8), allocatable :: actual(:), expected(:)
        integer(int32) :: i_byte

        call read_file_bytes(actual_path, actual)
        call read_file_bytes(expected_path, expected)
        call assert_equal_int(size(actual, kind=int32), size(expected, kind=int32), msg//": file size in bytes")
        if (size(actual) /= size(expected)) return
        do i_byte = 1, size(expected, kind=int32)
            if (actual(i_byte) /= expected(i_byte)) then
                call assert_equal_int(int(actual(i_byte), int32), int(expected(i_byte), int32), msg//": byte "//i_byte)
                return
            end if
        end do
    end subroutine check_same_file

    !> Writes the `i_bad`-th broken file of BAD_FILE_NAMES to `broken_path`, made from the valid
    !| file at `valid_path`, whose header has `n_extents` extents, the first `first_extent`;
    !| `other_type_code` is a type code the loader does not take.
    subroutine write_broken_file(i_bad, broken_path, valid_path, other_type_code, n_extents, first_extent)
        integer(int32), intent(in) :: i_bad, other_type_code, n_extents, first_extent
        character(*), intent(in) :: broken_path, valid_path
        integer(int8), allocatable :: valid(:)
        integer(int32) :: unit, status, header_size

        call read_file_bytes(valid_path, valid)
        header_size = 12 + 4*n_extents
        call open_for_writing(broken_path, unit, status)
        if (status /= 0) return
        select case (i_bad)
        case (1)
            ! nothing
        case (2)
            write (unit) valid(1:3)
        case (3)
            write (unit) MAGIC + 1_int32, valid(5:)
        case (4)
            write (unit) valid(1:4)
        case (5)
            write (unit) valid(1:8)
        case (6)
            write (unit) valid(1:8), 16_int32, valid(13:)
        case (7)
            write (unit) valid(1:8), -1_int32, valid(13:)
        case (8)
            write (unit) valid(1:header_size - 4)
        case (9)
            write (unit) valid(1:4), other_type_code, valid(9:)
        case (10)
            write (unit) valid(1:8), n_extents + 1_int32, valid(13:header_size), 1_int32, valid(header_size + 1:)
        case (11)
            write (unit) valid(1:12), first_extent + 1_int32, valid(17:), valid(header_size + 1:)
        case default
            write (unit) valid(1:size(valid) - 1)
        end select
        close (unit)
    end subroutine write_broken_file

    ! ---- the data, one set per kind ----

    !> Three ids in strings of 6: gene22 fills its string, the others are padded with blanks.
    subroutine three_strings(strings)
        character(len=6), intent(out) :: strings(3)

        strings = ["g1    ", "gene22", "x     "]
    end subroutine three_strings

    !> Four family indices, the sentinel 0 and the largest and smallest of the symmetric int32
    !| range among them.
    subroutine four_integers(integers)
        integer(int32), intent(out) :: integers(4)

        integers = [7_int32, 0_int32, huge(1_int32), -huge(1_int32)]
    end subroutine four_integers

    !> A 2 x 3 matrix whose values are exact in binary, the ends of the finite range and a
    !| negative zero among them.
    subroutine two_by_three_reals(reals)
        real(real64), intent(out) :: reals(2, 3)

        reals = reshape([1.5_real64, -2.0_real64, 0.25_real64, huge(1.0_real64), -tiny(1.0_real64), -0.0_real64], [2, 3])
    end subroutine two_by_three_reals

    !> NaN, +Inf and -Inf in a 1 x 3 matrix.
    subroutine non_finite_reals(reals)
        real(real64), intent(out) :: reals(1, 3)

        reals(1, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        reals(1, 2) = ieee_value(1.0_real64, ieee_positive_inf)
        reals(1, 3) = ieee_value(1.0_real64, ieee_negative_inf)
    end subroutine non_finite_reals

    ! ---- strings: gene_ids, family_ids ----

    !> The layout of three strings and of none, written by `save`.
    subroutine check_saves_strings(save, case_name)
        procedure(save_strings_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), no_strings(0)
        character(len=:), allocatable :: saved, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call three_strings(strings)
        call save(strings, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": three strings, ierr")
        call write_strings_file(expected, strings)
        call check_same_file(saved, expected, case_name//": three strings")

        call save(no_strings, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no string, ierr")
        call write_strings_file(expected, no_strings)
        call check_same_file(saved, expected, case_name//": no string")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_strings

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_strings(save, case_name)
        procedure(save_strings_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3)
        character(len=:), allocatable :: saved, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call write_scratch_file(saved, repeat("x", 100))
        call three_strings(strings)
        call save(strings, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call write_strings_file(expected, strings)
        call check_same_file(saved, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_strings

    !> Three strings and none, read by `load` from files in the layout.
    subroutine check_loads_strings(load, case_name)
        procedure(load_strings_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), loaded(3), no_strings(0)
        character(len=:), allocatable :: input
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        input = scratch_file(scratch, "input.bin")
        call three_strings(strings)
        call write_strings_file(input, strings)
        loaded = "unset"
        call load(loaded, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": three strings, ierr")
        call assert_equal_array_char(loaded, strings, 6, 3, case_name//": three strings")

        call write_strings_file(input, no_strings)
        call load(no_strings, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no string")
        call remove_scratch_directory(scratch)
    end subroutine check_loads_strings

    !> Every broken file of BAD_FILE_NAMES is refused by `load` with its code; a string length
    !| other than the array's (5 for strings of 6) counts as another type.
    subroutine check_load_refuses_strings(load, case_name)
        procedure(load_strings_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), loaded(3)
        character(len=:), allocatable :: valid, broken
        integer(int32) :: ierr, i_bad

        call create_scratch_directory(scratch, case_name)
        valid = scratch_file(scratch, "valid.bin")
        broken = scratch_file(scratch, "broken.bin")
        call three_strings(strings)
        call write_strings_file(valid, strings)
        do i_bad = 1, N_BAD_FILES
            call write_broken_file(i_bad, broken, valid, 5_int32, 1_int32, 3_int32)
            call load(loaded, broken, ierr)
            call assert_err(ierr, BAD_FILE_CODES(i_bad), case_name//": "//trim(BAD_FILE_NAMES(i_bad)))
        end do
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_strings

    !> What `save` writes, `load` reads back unchanged.
    subroutine check_round_trip_strings(save, load, case_name)
        procedure(save_strings_interface) :: save
        procedure(load_strings_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), loaded(3)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call three_strings(strings)
        call save(strings, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": save")
        loaded = "unset"
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": load")
        call assert_equal_array_char(loaded, strings, 6, 3, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_round_trip_strings

    ! ---- integers: gene_to_family ----

    !> The layout of four integers and of none, written by `save`.
    subroutine check_saves_integers(save, case_name)
        procedure(save_integers_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), no_integers(0), ierr
        character(len=:), allocatable :: saved, expected

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call four_integers(integers)
        call save(integers, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": four integers, ierr")
        call write_integers_file(expected, integers)
        call check_same_file(saved, expected, case_name//": four integers")

        call save(no_integers, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no integer, ierr")
        call write_integers_file(expected, no_integers)
        call check_same_file(saved, expected, case_name//": no integer")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_integers

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_integers(save, case_name)
        procedure(save_integers_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), ierr
        character(len=:), allocatable :: saved, expected

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call write_scratch_file(saved, repeat("x", 100))
        call four_integers(integers)
        call save(integers, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call write_integers_file(expected, integers)
        call check_same_file(saved, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_integers

    !> Four integers and none, read by `load` from files in the layout.
    subroutine check_loads_integers(load, case_name)
        procedure(load_integers_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), loaded(4), no_integers(0), ierr
        character(len=:), allocatable :: input

        call create_scratch_directory(scratch, case_name)
        input = scratch_file(scratch, "input.bin")
        call four_integers(integers)
        call write_integers_file(input, integers)
        loaded = -5
        call load(loaded, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": four integers, ierr")
        call assert_equal_array_int(loaded, integers, 4, case_name//": four integers")

        call write_integers_file(input, no_integers)
        call load(no_integers, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no integer")
        call remove_scratch_directory(scratch)
    end subroutine check_loads_integers

    !> Every broken file of BAD_FILE_NAMES is refused by `load` with its code; a real file is
    !| another type.
    subroutine check_load_refuses_integers(load, case_name)
        procedure(load_integers_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), loaded(4), ierr, i_bad
        character(len=:), allocatable :: valid, broken

        call create_scratch_directory(scratch, case_name)
        valid = scratch_file(scratch, "valid.bin")
        broken = scratch_file(scratch, "broken.bin")
        call four_integers(integers)
        call write_integers_file(valid, integers)
        do i_bad = 1, N_BAD_FILES
            call write_broken_file(i_bad, broken, valid, REAL_CODE, 1_int32, 4_int32)
            call load(loaded, broken, ierr)
            call assert_err(ierr, BAD_FILE_CODES(i_bad), case_name//": "//trim(BAD_FILE_NAMES(i_bad)))
        end do
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_integers

    !> What `save` writes, `load` reads back unchanged.
    subroutine check_round_trip_integers(save, load, case_name)
        procedure(save_integers_interface) :: save
        procedure(load_integers_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), loaded(4), ierr
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call four_integers(integers)
        call save(integers, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": save")
        loaded = -5
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": load")
        call assert_equal_array_int(loaded, integers, 4, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_round_trip_integers

    ! ---- reals: expression_vectors, family_centroids, shift_vectors ----
    ! Loaded reals are compared as bytes, by writing them in the layout and comparing that file
    ! with the one they came from, so the sign of a zero and the bits of a NaN count.

    !> The layout of a 2 x 3 matrix and of a 0 x 3 one, written by `save`.
    subroutine check_saves_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), no_rows(0, 3)
        character(len=:), allocatable :: saved, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call two_by_three_reals(reals)
        call save(reals, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 2 x 3, ierr")
        call write_reals_file(expected, reals)
        call check_same_file(saved, expected, case_name//": 2 x 3")

        call save(no_rows, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 0 x 3, ierr")
        call write_reals_file(expected, no_rows)
        call check_same_file(saved, expected, case_name//": 0 x 3")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_reals

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3)
        character(len=:), allocatable :: saved, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call write_scratch_file(saved, repeat("x", 100))
        call two_by_three_reals(reals)
        call save(reals, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call write_reals_file(expected, reals)
        call check_same_file(saved, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_reals

    !> A 2 x 3 matrix and a 0 x 3 one, read by `load` from files in the layout.
    subroutine check_loads_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), loaded(2, 3), no_rows(0, 3)
        character(len=:), allocatable :: input, read_back
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        input = scratch_file(scratch, "input.bin")
        read_back = scratch_file(scratch, "read_back.bin")
        call two_by_three_reals(reals)
        call write_reals_file(input, reals)
        loaded = -5.0_real64
        call load(loaded, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 2 x 3, ierr")
        call write_reals_file(read_back, loaded)
        call check_same_file(read_back, input, case_name//": 2 x 3")

        call write_reals_file(input, no_rows)
        call load(no_rows, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 0 x 3")
        call remove_scratch_directory(scratch)
    end subroutine check_loads_reals

    !> Every broken file of BAD_FILE_NAMES is refused by `load` with its code; an integer file is
    !| another type.
    subroutine check_load_refuses_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), loaded(2, 3)
        character(len=:), allocatable :: valid, broken
        integer(int32) :: ierr, i_bad

        call create_scratch_directory(scratch, case_name)
        valid = scratch_file(scratch, "valid.bin")
        broken = scratch_file(scratch, "broken.bin")
        call two_by_three_reals(reals)
        call write_reals_file(valid, reals)
        do i_bad = 1, N_BAD_FILES
            call write_broken_file(i_bad, broken, valid, INTEGER_CODE, 2_int32, 2_int32)
            call load(loaded, broken, ierr)
            call assert_err(ierr, BAD_FILE_CODES(i_bad), case_name//": "//trim(BAD_FILE_NAMES(i_bad)))
        end do
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_reals

    !> What `save` writes, `load` reads back unchanged, bit for bit.
    subroutine check_round_trip_reals(save, load, case_name)
        procedure(save_reals_interface) :: save
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), loaded(2, 3)
        character(len=:), allocatable :: saved, read_back, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        read_back = scratch_file(scratch, "read_back.bin")
        expected = scratch_file(scratch, "expected.bin")
        call two_by_three_reals(reals)
        call save(reals, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": save")
        loaded = -5.0_real64
        call load(loaded, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": load")
        call write_reals_file(read_back, loaded)
        call write_reals_file(expected, reals)
        call check_same_file(read_back, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_round_trip_reals

    !> `save` writes NaN, +Inf and -Inf bit for bit, without an error, as the module documents.
    subroutine check_saves_non_finite_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(1, 3)
        character(len=:), allocatable :: saved, expected
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        saved = scratch_file(scratch, "saved.bin")
        expected = scratch_file(scratch, "expected.bin")
        call non_finite_reals(reals)
        call save(reals, saved, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call write_reals_file(expected, reals)
        call check_same_file(saved, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_saves_non_finite_reals

    !> `load` reads NaN, +Inf and -Inf bit for bit, without an error, as the module documents.
    subroutine check_loads_non_finite_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(1, 3), loaded(1, 3)
        character(len=:), allocatable :: input, read_back
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        input = scratch_file(scratch, "input.bin")
        read_back = scratch_file(scratch, "read_back.bin")
        call non_finite_reals(reals)
        call write_reals_file(input, reals)
        call load(loaded, input, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call write_reals_file(read_back, loaded)
        call check_same_file(read_back, input, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_loads_non_finite_reals
    !> A path in a directory that does not exist is refused by `save` with ERR_FILE_OPEN, and the
    !| error unit stays connected.
    subroutine check_save_refuses_an_unwritable_path_strings(save, case_name)
        procedure(save_strings_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: values(3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call three_strings(values)
        call clear_stack()
        call save(values, scratch%path//"/no_such_directory/file.bin", ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_refuses_an_unwritable_path_strings

    !> A file that does not exist is refused by `load` with ERR_FILE_OPEN, and the error unit stays
    !| connected.
    subroutine check_load_refuses_a_missing_file_strings(load, case_name)
        procedure(load_strings_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: values(3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call clear_stack()
        call load(values, scratch_file(scratch, "missing.bin"), ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_a_missing_file_strings

    !> A path in a directory that does not exist is refused by `save` with ERR_FILE_OPEN, and the
    !| error unit stays connected.
    subroutine check_save_refuses_an_unwritable_path_integers(save, case_name)
        procedure(save_integers_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: values(4)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call four_integers(values)
        call clear_stack()
        call save(values, scratch%path//"/no_such_directory/file.bin", ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_refuses_an_unwritable_path_integers

    !> A file that does not exist is refused by `load` with ERR_FILE_OPEN, and the error unit stays
    !| connected.
    subroutine check_load_refuses_a_missing_file_integers(load, case_name)
        procedure(load_integers_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: values(4)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call clear_stack()
        call load(values, scratch_file(scratch, "missing.bin"), ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_a_missing_file_integers

    !> A path in a directory that does not exist is refused by `save` with ERR_FILE_OPEN, and the
    !| error unit stays connected.
    subroutine check_save_refuses_an_unwritable_path_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: values(2, 3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call two_by_three_reals(values)
        call clear_stack()
        call save(values, scratch%path//"/no_such_directory/file.bin", ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_refuses_an_unwritable_path_reals

    !> A file that does not exist is refused by `load` with ERR_FILE_OPEN, and the error unit stays
    !| connected.
    subroutine check_load_refuses_a_missing_file_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: values(2, 3)
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        call clear_stack()
        call load(values, scratch_file(scratch, "missing.bin"), ierr)
        call assert_err(ierr, ERR_FILE_OPEN, case_name)
        call check_error_unit_connected(case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_load_refuses_a_missing_file_reals
end module mod_test_tox_data_read_write_fixtures
