!> What the children of the `tox_data_read_write` suite share: the published file layout, written
!| out by hand, and the checks every `save_*` and `load_*` procedure gets, once per kind of array.
!|
!| A file is an unformatted stream: the magic number 1178677808 (`'FA20'`, hex 46413230), the
!| type code (-1 integer, -2 real, a string length for strings), the number of dimensions, the
!| extents, then the elements in column-major order (f42_serde_arrays_utils documents it). Each
!| `save_*` case compares a written file byte for byte with that layout, and each `load_*` case
!| reads a file written byte for byte by the case, so neither direction is tested only through
!| the other.
!|
!| The procedures of one kind take the same arguments, so each check takes the procedure under
!| test as an argument: `tox_data_read_write`'s three real matrices share one set of checks, its
!| two string arrays another. Each check is one case of a child, in a scratch directory of its
!| own, named after the case.
!|
!| Not checked here: a missing file, or a path that cannot be written. Both make the serde layer
!| close a unit number it never opened, which closes stderr or aborts the process (see the batch's
!| report); they cannot be cases until that is fixed.
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

    !> The magic number every array file starts with, 'FA20' read as an int32.
    integer(int32), parameter :: MAGIC = 1178677808
    integer(int32), parameter :: INTEGER_CODE = -1, REAL_CODE = -2
    !> The MOLD of a `transfer` to bytes.
    integer(int8), parameter :: BYTE_MOLD(1) = 0_int8

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

    ! ---- the layout ----

    !> The bytes of a header: magic number, type code, number of dimensions, extents.
    function header_bytes(type_code, extents) result(bytes)
        integer(int32), intent(in) :: type_code, extents(:)
        integer(int8), allocatable :: bytes(:)
        integer(int32) :: words(3 + size(extents))

        words(1:3) = [MAGIC, type_code, size(extents)]
        words(4:) = extents
        bytes = transfer(words, BYTE_MOLD)
    end function header_bytes

    !> The bytes of `file_path`, all of them; none for a file that cannot be read.
    function file_bytes(file_path) result(bytes)
        character(*), intent(in) :: file_path
        integer(int8), allocatable :: bytes(:)
        integer(int32) :: unit, status, file_size

        allocate (bytes(0))
        open (newunit=unit, file=file_path, access='stream', form='unformatted', status='old', action='read', &
              iostat=status)
        if (status /= 0) return
        inquire (unit=unit, size=file_size)
        deallocate (bytes)
        allocate (bytes(file_size))
        if (file_size > 0) read (unit, iostat=status) bytes
        close (unit)
    end function file_bytes

    !> Writes exactly `bytes` to `file_path`.
    subroutine write_bytes(file_path, bytes)
        character(*), intent(in) :: file_path
        integer(int8), intent(in) :: bytes(:)
        integer(int32) :: unit, status

        open (newunit=unit, file=file_path, access='stream', form='unformatted', status='replace', action='write', &
              iostat=status)
        call assert_equal_int(status, 0_int32, "write_bytes: open "//file_path)
        if (status /= 0) return
        if (size(bytes) > 0) write (unit) bytes
        close (unit)
    end subroutine write_bytes

    !> 100 bytes, more than any file a case saves, to be replaced by one.
    function longer_file() result(bytes)
        integer(int8), allocatable :: bytes(:)

        bytes = transfer(repeat("x", 100), BYTE_MOLD)
    end function longer_file

    !> Fails unless `file_path` holds exactly `expected`.
    subroutine check_file_holds(file_path, expected, msg)
        character(*), intent(in) :: file_path, msg
        integer(int8), intent(in) :: expected(:)
        integer(int8), allocatable :: actual(:)
        integer(int32), allocatable :: actual_words(:), expected_words(:)

        actual = file_bytes(file_path)
        call assert_equal_int(size(actual), size(expected), msg//": file size in bytes")
        if (size(actual) /= size(expected)) return
        actual_words = int(actual, int32)
        expected_words = int(expected, int32)
        call assert_equal_array_int(actual_words, expected_words, size(expected), msg//": bytes")
    end subroutine check_file_holds

    !> The `i_bad`-th broken file of BAD_FILE_NAMES, made from a valid file's `type_code`,
    !| `extents` and `payload`; `other_type_code` is a type code the loader does not take.
    function bad_file(i_bad, type_code, other_type_code, extents, payload) result(bytes)
        integer(int32), intent(in) :: i_bad, type_code, other_type_code, extents(:)
        integer(int8), intent(in) :: payload(:)
        integer(int8), allocatable :: bytes(:), valid_header(:)
        integer(int32) :: wrong_extents(size(extents)), extra_rank(size(extents) + 1)

        valid_header = header_bytes(type_code, extents)
        select case (i_bad)
        case (1)
            allocate (bytes(0))
        case (2)
            bytes = valid_header(1:3)
        case (3)
            bytes = [transfer(MAGIC + 1_int32, BYTE_MOLD), valid_header(5:), payload]
        case (4)
            bytes = valid_header(1:4)
        case (5)
            bytes = valid_header(1:8)
        case (6)
            bytes = [valid_header(1:8), transfer(16_int32, BYTE_MOLD), valid_header(13:), payload]
        case (7)
            bytes = [valid_header(1:8), transfer(-1_int32, BYTE_MOLD), valid_header(13:), payload]
        case (8)
            bytes = valid_header(1:size(valid_header) - 4)
        case (9)
            bytes = [header_bytes(other_type_code, extents), payload]
        case (10)
            extra_rank(1:size(extents)) = extents
            extra_rank(size(extents) + 1) = 1
            bytes = [header_bytes(type_code, extra_rank), payload]
        case (11)
            wrong_extents = extents
            wrong_extents(1) = wrong_extents(1) + 1
            bytes = [header_bytes(type_code, wrong_extents), payload, payload]
        case default
            bytes = [valid_header, payload(1:size(payload) - 1)]
        end select
    end function bad_file

    ! ---- the data, one set per kind ----

    !> Three ids in strings of 6: gene22 fills its string, the others are padded with blanks.
    subroutine three_strings(strings, payload)
        character(len=6), intent(out) :: strings(3)
        integer(int8), allocatable, intent(out) :: payload(:)
        character(len=18) :: written

        strings = ["g1    ", "gene22", "x     "]
        written = "g1    gene22x     "
        payload = transfer(written, BYTE_MOLD)
    end subroutine three_strings

    !> Four family indices, the sentinel 0 and the largest and smallest of the symmetric int32
    !| range among them.
    subroutine four_integers(integers, payload)
        integer(int32), intent(out) :: integers(4)
        integer(int8), allocatable, intent(out) :: payload(:)

        integers = [7_int32, 0_int32, huge(1_int32), -huge(1_int32)]
        payload = transfer(integers, BYTE_MOLD)
    end subroutine four_integers

    !> A 2 x 3 matrix whose values are exact in binary, the ends of the finite range and a
    !| negative zero among them.
    subroutine two_by_three_reals(reals, payload)
        real(real64), intent(out) :: reals(2, 3)
        integer(int8), allocatable, intent(out) :: payload(:)

        reals = reshape([1.5_real64, -2.0_real64, 0.25_real64, huge(1.0_real64), -tiny(1.0_real64), -0.0_real64], [2, 3])
        payload = transfer(reals, BYTE_MOLD)
    end subroutine two_by_three_reals

    !> NaN, +Inf and -Inf in a 1 x 3 matrix.
    subroutine non_finite_reals(reals, payload)
        real(real64), intent(out) :: reals(1, 3)
        integer(int8), allocatable, intent(out) :: payload(:)

        reals(1, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        reals(1, 2) = ieee_value(1.0_real64, ieee_positive_inf)
        reals(1, 3) = ieee_value(1.0_real64, ieee_negative_inf)
        payload = transfer(reals, BYTE_MOLD)
    end subroutine non_finite_reals

    !> Fails unless the bytes of `reals` are `expected`: compared as bytes, so the sign of a zero
    !| and the bits of a NaN count.
    subroutine check_real_bytes(reals, expected, msg)
        real(real64), intent(in) :: reals(:, :)
        integer(int8), intent(in) :: expected(:)
        character(*), intent(in) :: msg
        integer(int32), allocatable :: actual_words(:), expected_words(:)

        actual_words = int(transfer(reals, BYTE_MOLD), int32)
        expected_words = int(expected, int32)
        call assert_equal_int(size(actual_words), size(expected_words), msg//": size in bytes")
        if (size(actual_words) /= size(expected_words)) return
        call assert_equal_array_int(actual_words, expected_words, size(expected_words), msg)
    end subroutine check_real_bytes

    ! ---- strings: gene_ids, family_ids ----

    !> The layout of three strings and of none, written by `save`.
    subroutine check_saves_strings(save, case_name)
        procedure(save_strings_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), no_strings(0)
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call three_strings(strings, payload)
        call save(strings, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": three strings, ierr")
        expected = [header_bytes(6_int32, [3_int32]), payload]
        call check_file_holds(file_path, expected, case_name//": three strings")

        call save(no_strings, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no string, ierr")
        expected = header_bytes(6_int32, [0_int32])
        call check_file_holds(file_path, expected, case_name//": no string")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_strings

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_strings(save, case_name)
        procedure(save_strings_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3)
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call write_bytes(file_path, longer_file())
        call three_strings(strings, payload)
        call save(strings, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        expected = [header_bytes(6_int32, [3_int32]), payload]
        call check_file_holds(file_path, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_strings

    !> Three strings and none, read by `load` from files written byte for byte.
    subroutine check_loads_strings(load, case_name)
        procedure(load_strings_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        character(len=6) :: strings(3), loaded(3), no_strings(0)
        integer(int8), allocatable :: payload(:), written(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call three_strings(strings, payload)
        written = [header_bytes(6_int32, [3_int32]), payload]
        call write_bytes(file_path, written)
        loaded = "unset"
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": three strings, ierr")
        call assert_equal_array_char(loaded, strings, 6, 3, case_name//": three strings")

        call write_bytes(file_path, header_bytes(6_int32, [0_int32]))
        call load(no_strings, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr, i_bad

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call three_strings(strings, payload)
        do i_bad = 1, N_BAD_FILES
            call write_bytes(file_path, bad_file(i_bad, 6_int32, 5_int32, [3_int32], payload))
            call load(loaded, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "strings.bin")
        call three_strings(strings, payload)
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
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call four_integers(integers, payload)
        call save(integers, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": four integers, ierr")
        expected = [header_bytes(INTEGER_CODE, [4_int32]), payload]
        call check_file_holds(file_path, expected, case_name//": four integers")

        call save(no_integers, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": no integer, ierr")
        expected = header_bytes(INTEGER_CODE, [0_int32])
        call check_file_holds(file_path, expected, case_name//": no integer")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_integers

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_integers(save, case_name)
        procedure(save_integers_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), ierr
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call write_bytes(file_path, longer_file())
        call four_integers(integers, payload)
        call save(integers, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        expected = [header_bytes(INTEGER_CODE, [4_int32]), payload]
        call check_file_holds(file_path, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_integers

    !> Four integers and none, read by `load` from files written byte for byte.
    subroutine check_loads_integers(load, case_name)
        procedure(load_integers_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        integer(int32) :: integers(4), loaded(4), no_integers(0), ierr
        integer(int8), allocatable :: payload(:), written(:)
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call four_integers(integers, payload)
        written = [header_bytes(INTEGER_CODE, [4_int32]), payload]
        call write_bytes(file_path, written)
        loaded = -5
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": four integers, ierr")
        call assert_equal_array_int(loaded, integers, 4, case_name//": four integers")

        call write_bytes(file_path, header_bytes(INTEGER_CODE, [0_int32]))
        call load(no_integers, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call four_integers(integers, payload)
        do i_bad = 1, N_BAD_FILES
            call write_bytes(file_path, bad_file(i_bad, INTEGER_CODE, REAL_CODE, [4_int32], payload))
            call load(loaded, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "integers.bin")
        call four_integers(integers, payload)
        call save(integers, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": save")
        loaded = -5
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": load")
        call assert_equal_array_int(loaded, integers, 4, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_round_trip_integers

    ! ---- reals: expression_vectors, family_centroids, shift_vectors ----

    !> The layout of a 2 x 3 matrix and of a 0 x 3 one, written by `save`.
    subroutine check_saves_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), no_rows(0, 3)
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call two_by_three_reals(reals, payload)
        call save(reals, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 2 x 3, ierr")
        expected = [header_bytes(REAL_CODE, [2_int32, 3_int32]), payload]
        call check_file_holds(file_path, expected, case_name//": 2 x 3")

        call save(no_rows, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 0 x 3, ierr")
        expected = header_bytes(REAL_CODE, [0_int32, 3_int32])
        call check_file_holds(file_path, expected, case_name//": 0 x 3")
        call remove_scratch_directory(scratch)
    end subroutine check_saves_reals

    !> A longer file in the way is replaced by what `save` writes, not overwritten in part.
    subroutine check_save_replaces_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3)
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call write_bytes(file_path, longer_file())
        call two_by_three_reals(reals, payload)
        call save(reals, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        expected = [header_bytes(REAL_CODE, [2_int32, 3_int32]), payload]
        call check_file_holds(file_path, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_save_replaces_reals

    !> A 2 x 3 matrix and a 0 x 3 one, read by `load` from files written byte for byte.
    subroutine check_loads_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(2, 3), loaded(2, 3), no_rows(0, 3)
        integer(int8), allocatable :: payload(:), written(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call two_by_three_reals(reals, payload)
        written = [header_bytes(REAL_CODE, [2_int32, 3_int32]), payload]
        call write_bytes(file_path, written)
        loaded = -5.0_real64
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": 2 x 3, ierr")
        call check_real_bytes(loaded, payload, case_name//": 2 x 3")

        call write_bytes(file_path, header_bytes(REAL_CODE, [0_int32, 3_int32]))
        call load(no_rows, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr, i_bad

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call two_by_three_reals(reals, payload)
        do i_bad = 1, N_BAD_FILES
            call write_bytes(file_path, bad_file(i_bad, REAL_CODE, INTEGER_CODE, [2_int32, 3_int32], payload))
            call load(loaded, file_path, ierr)
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
        integer(int8), allocatable :: payload(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call two_by_three_reals(reals, payload)
        call save(reals, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": save")
        loaded = -5.0_real64
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": load")
        call check_real_bytes(loaded, payload, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_round_trip_reals

    !> Pins today's behaviour: `save` writes NaN, +Inf and -Inf bit for bit, without an error.
    !| The procedures do not document whether a stored matrix may hold them (the batch's report
    !| asks).
    subroutine check_saves_non_finite_reals(save, case_name)
        procedure(save_reals_interface) :: save
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(1, 3)
        integer(int8), allocatable :: payload(:), expected(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call non_finite_reals(reals, payload)
        call save(reals, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        expected = [header_bytes(REAL_CODE, [1_int32, 3_int32]), payload]
        call check_file_holds(file_path, expected, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_saves_non_finite_reals

    !> Pins today's behaviour: `load` reads NaN, +Inf and -Inf bit for bit, without an error (see
    !| `check_saves_non_finite_reals`).
    subroutine check_loads_non_finite_reals(load, case_name)
        procedure(load_reals_interface) :: load
        character(*), intent(in) :: case_name
        type(scratch_directory_type) :: scratch
        real(real64) :: reals(1, 3), loaded(1, 3)
        integer(int8), allocatable :: payload(:), written(:)
        character(len=:), allocatable :: file_path
        integer(int32) :: ierr

        call create_scratch_directory(scratch, case_name)
        file_path = scratch_file(scratch, "reals.bin")
        call non_finite_reals(reals, payload)
        written = [header_bytes(REAL_CODE, [1_int32, 3_int32]), payload]
        call write_bytes(file_path, written)
        call load(loaded, file_path, ierr)
        call assert_err(ierr, ERR_OK, case_name//": ierr")
        call check_real_bytes(loaded, payload, case_name)
        call remove_scratch_directory(scratch)
    end subroutine check_loads_non_finite_reals
end module mod_test_tox_data_read_write_fixtures
