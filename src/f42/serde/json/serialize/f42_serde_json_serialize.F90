#include <src/macros.h>

!> AUTHOR_FRANZ_ERIC_SILL
!| A streaming JSON writer: a document is written to its file while it is being described, so
!| writing it needs no tree, no copy of the data and no memory beyond one internal buffer.
!|
!| ### Streaming model
!|
!| A [[f42_serde_json_serialize(module):json_writer(type)]] is opened on a file with
!| [[f42_serde_json_serialize(module):json_open(subroutine)]], receives the document as a
!| sequence of calls, and is finished with
!| [[f42_serde_json_serialize(module):json_close(subroutine)]]:
!|
!| - [[f42_serde_json_serialize(module):json_begin_object(subroutine)]] and
!|   [[f42_serde_json_serialize(module):json_end_object(subroutine)]] open and close an object,
!|   [[f42_serde_json_serialize(module):json_begin_array(subroutine)]] and
!|   [[f42_serde_json_serialize(module):json_end_array(subroutine)]] an array;
!| - [[f42_serde_json_serialize(module):json_member(interface)]] writes one `"key":value` pair
!|   inside an object, [[f42_serde_json_serialize(module):json_element(interface)]] one value
!|   inside an array (or as the root), and
!|   [[f42_serde_json_serialize(module):json_null(subroutine)]] writes `null` in either place.
!|
!| Both generics take their value as `class(*)`, a scalar or a rank-1 array (written as one
!| compact JSON array), so one call serves every type. The supported types are
!| `integer(int32)`, `integer(int64)`, `real(real64)`, `logical(c_bool)`, the default `logical`
!| and `character(*)`. Any other type -- `real(real32)`, `integer(int8)`, `complex`, a derived
!| type -- is refused at run time with `ERR_INVALID_INPUT`, never written as `null`. An empty
!| array is written as `[]` whatever its element type, so an empty array of an otherwise
!| unsupported type is accepted too (there is nothing of that type to write). nvfortran
!| needs version 26.9 or later: 25.7 passes a character expression such as `text//"x"` to a
!| `class(*)` argument as an empty string, and crashes in its copy-in at `-O3`.
!|
!| A container opened with a key is a member of the enclosing object; one opened without a key
!| is an element of the enclosing array, or the root. The whole file is one line without
!| whitespace, ended by a single line feed.
!|
!| ### What the writer refuses
!|
!| A nesting state machine checks every call, so a call sequence can only ever produce valid
!| JSON. It refuses, with `ERR_JSON_SEQUENCE` from [[tox_errors(module)]]:
!|
!| - a value without a key inside an object, and a value with a key inside an array or as the
!|   root;
!| - an end that does not match the innermost open container, or one with nothing open;
!| - anything after the root value is complete: a document has exactly one root, and any
!|   value may be it;
!| - closing the writer while a container is still open or before a root was written;
!| - any call on a writer that is not open, and opening a writer that already is.
!|
!| Duplicate keys in one object are **not** detected: that would mean remembering every key,
!| which the streaming model exists to avoid. Writing each key once is the caller's
!| responsibility. Nesting has no depth limit.
!|
!| ### Errors are sticky
!|
!| Only `json_open` and `json_close` return an error code; every other call reports through
!| the writer. The first failure is kept in `w%ierr`, and `w%failing_call` numbers the call that
!| caused it, counting `json_open` as call 1. Every later call is a no-op, so a sequence of
!| calls can be written without checking each one and `json_close` returns the first failure.
!| On **any** failure the file is deleted at once: a file the writer leaves behind is always a
!| complete, valid document. The codes are `ERR_FILE_OPEN` (the file could not be created,
!| including when it already exists), `ERR_WRITE_DATA`, `ERR_JSON_SEQUENCE`, `ERR_NAN_INF`,
!| `ERR_INVALID_UTF8`, `ERR_INVALID_INPUT` (a value of an unsupported type) and `ERR_ALLOC_FAIL`.
!| A call is checked against the nesting rules before its value is looked at.
!|
!| ### Numbers
!|
!| Integers are written in full (`I0`). A `real(real64)` is written with 17 significant digits
!| (`ES24.16E3`, leading blank removed, as in `-1.2345678901234567E-005`), which reads back as
!| exactly the same double. A zero is recognised by its IEEE class, never by comparing it with
!| 0, and written as `0.0` or `-0.0`, so a negative zero keeps its sign and a subnormal number
!| is never mistaken for zero, not even where denormals are flushed. NaN and infinity have no
!| JSON form and are refused with `ERR_NAN_INF`.
!|
!| ### Strings
!|
!| Keys and string values follow the same rules. Trailing blanks are removed, since Fortran
!| strings are blank-padded; leading blanks and every other byte are kept, so a value that
!| genuinely ends in a blank cannot be written. `"` and `\` are escaped, as are the control
!| characters: `\b`, `\f`, `\n`, `\r` and `\t` in their short forms, every other byte below
!| 0x20 (NUL included) as `\u00xx` with lowercase hex digits. `/` and DEL are written as
!| they are. Bytes from 0x80 up must form valid UTF-8 and are copied unchanged. Anything else
!| (a stray continuation byte, an overlong form, a surrogate U+D800 to U+DFFF, a code point
!| above U+10FFFF or a truncated sequence) is refused with `ERR_INVALID_UTF8`, because a JSON
!| file must be UTF-8. [[f42_serde_json_serialize(module):json_is_valid_utf8(function)]]
!| applies the same test, so a caller can validate its strings before it opens a file.
!|
!| ### Input/output
!|
!| The file is written as an unformatted stream through a 256 KiB buffer the writer owns,
!| which is the one mode that is byte-exact on every supported compiler. Every write and the
!| close are checked, and after closing, the file's size on disk is compared with the number
!| of bytes written, because a runtime may lose an error such as a full disk on a small write.
module f42_serde_json_serialize
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_class, ieee_class_type, operator(==), &
                                             ieee_positive_zero, ieee_negative_zero
    use tox_errors, only: set_ok, set_err, is_ok, is_err, &
                          ERR_OK, ERR_FILE_OPEN, ERR_WRITE_DATA, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          ERR_JSON_SEQUENCE, ERR_INVALID_UTF8, ERR_INVALID_INPUT
    M_IMPLICIT_NONE

    private
    public :: json_writer
    public :: json_open, json_close
    public :: json_begin_object, json_end_object, json_begin_array, json_end_array
    public :: json_member, json_element, json_null
    public :: json_is_valid_utf8

    integer(int8), parameter :: CONTAINER_OBJECT = 1_int8
        !! Marks an open object on the container stack
    integer(int8), parameter :: CONTAINER_ARRAY = 2_int8
        !! Marks an open array on the container stack
    character(len=*), parameter :: OPENING_BRACKETS = '{['
        !! The opening bracket of each container kind, indexed by `CONTAINER_OBJECT` or `CONTAINER_ARRAY`
    character(len=*), parameter :: CLOSING_BRACKETS = '}]'
        !! The closing bracket of each container kind, indexed like `OPENING_BRACKETS`
    integer(int32), parameter :: INITIAL_STACK_CAPACITY = 16
        !! Container stack slots allocated by `json_open`; the stack doubles when full
    integer(int32), parameter :: BUFFER_CAPACITY = 262144
        !! Size of the output buffer in bytes (256 KiB)
    integer(int32), parameter :: NUMBER_BLOCK = 256
        !! Numbers formatted per internal write when writing an array
    integer(int32), parameter :: REAL_FIELD_WIDTH = 24
        !! Width of one `ES24.16E3` field
    integer(int32), parameter :: INTEGER_FIELD_WIDTH = 21
        !! Room for one integer (at most 20 characters) and its comma
    character(len=*), parameter :: REAL_FORMAT = '(*(ES24.16E3))'
        !! Format of one real or a block of them, each in a field of `REAL_FIELD_WIDTH`
    character(len=*), parameter :: INTEGER_FORMAT = '(*(I0,:,","))'
        !! Format of one integer or a block of them, comma-separated
    character(len=*), parameter :: BACKSLASH = achar(92)
        !! The backslash, spelled so that no compiler reads it as an escape in a literal
    character(len=*), parameter :: HEX_DIGITS = '0123456789abcdef'
        !! Lowercase hexadecimal digits for `\u00xx`

    !> A JSON document being written to a file; see [[f42_serde_json_serialize(module)]].
    !|
    !| Declare one, pass it to `json_open`, describe the document, and pass it to `json_close`.
    !| A closed writer may be opened again.
    type :: json_writer
        private
        integer(int32), public :: ierr = ERR_OK
            !! The first error, `ERR_OK` while there is none. Later calls do not change it
        integer(int64), public :: failing_call = 0
            !! Number of the call that failed, counting `json_open` as call 1; 0 while there is no error
        logical(c_bool) :: is_open = .false.
            !! Whether the file is open (a `newunit` number is negative, so the unit cannot say)
        integer(int32) :: unit = 0
            !! Unit of the open file
        character(len=:), allocatable :: filename
            !! Name of the open file, kept to check its size after closing and to delete it
        character(len=:), allocatable :: buffer
            !! Output buffer. Allocatable on purpose: a fixed-size component makes a compiler place a local writer in static storage
        integer(int32) :: buffer_used = 0
            !! Bytes of `buffer` holding output not yet written
        integer(int64) :: bytes_written = 0
            !! Bytes handed to the file so far, compared with its size after closing
        integer(int64) :: call_count = 0
            !! Calls made since `json_open`, which is call 1
        integer(int8), allocatable :: container_kinds(:)
            !! The open containers, innermost last
        integer(int32) :: depth = 0
            !! Number of open containers
        logical(c_bool) :: container_is_empty = .true.
            !! Whether the innermost open container has no value yet, so the next one needs no comma
        logical(c_bool) :: root_written = .false.
            !! Whether the root value is complete, after which nothing may follow
    end type json_writer

    !> Write one `"key":value` pair inside an object.
    !|
    !| `call json_member(w, key, value)`: `key` comes first and is required. `value` is a scalar
    !| or a rank-1 array (written as a JSON array) of one of the supported types:
    !| `integer(int32)`, `integer(int64)`, `real(real64)`, `logical(c_bool)`, default `logical`
    !| and `character(*)`. Any other type fails the writer with `ERR_INVALID_INPUT`.
    interface json_member
        module procedure member_value, member_values
    end interface json_member

    !> Write one value inside an array, or as the root of the document.
    !|
    !| `call json_element(w, value)`, with `value` a scalar or a rank-1 array of a type
    !| [[f42_serde_json_serialize(module):json_member(interface)]] supports.
    interface json_element
        module procedure element_value, element_values
    end interface json_element


contains

    ! ============================================================================================
    ! Opening and closing
    ! ============================================================================================

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Create `filename` and open the writer on it.
    !|
    !| The file must not exist yet: an existing file is never overwritten, and is refused with
    !| `ERR_FILE_OPEN`, as is a path that cannot be created. Opening a writer that is already
    !| open fails that writer with `ERR_JSON_SEQUENCE` (and deletes its file).
    subroutine json_open(w, filename, ierr)
        type(json_writer), intent(inout) :: w
            !! The writer; it must not be open
        character(len=*), intent(in) :: filename
            !! Path of the file to create. Trailing blanks are ignored
        integer(int32), intent(out) :: ierr
            !! Error code: `ERR_OK`, `ERR_FILE_OPEN`, `ERR_ALLOC_FAIL` or `ERR_JSON_SEQUENCE`

        integer(int32) :: alloc_stat, iostat

        if (w%is_open) then
            w%call_count = w%call_count + 1
            call fail(w, ERR_JSON_SEQUENCE)
            ierr = w%ierr
            return
        end if

        call set_ok(w%ierr)
        w%failing_call = 0
        w%call_count = 1
        w%buffer_used = 0
        w%bytes_written = 0
        w%depth = 0
        w%container_is_empty = .true.
        w%root_written = .false.
        if (allocated(w%filename)) deallocate (w%filename)
        allocate (character(len=len_trim(filename)) :: w%filename, stat=alloc_stat)
        if (alloc_stat /= 0) then
            call fail(w, ERR_ALLOC_FAIL)
        else
            w%filename(:) = filename(1:len_trim(filename))
        end if

        if (.not. allocated(w%buffer) .and. is_ok(w%ierr)) then
            allocate (character(len=BUFFER_CAPACITY) :: w%buffer, stat=alloc_stat)
            if (alloc_stat /= 0) call fail(w, ERR_ALLOC_FAIL)
        end if
        if (.not. allocated(w%container_kinds) .and. is_ok(w%ierr)) then
            allocate (w%container_kinds(INITIAL_STACK_CAPACITY), stat=alloc_stat)
            if (alloc_stat /= 0) call fail(w, ERR_ALLOC_FAIL)
        end if

        if (is_ok(w%ierr)) then
            open (newunit=w%unit, file=w%filename, access='stream', form='unformatted', &
                  status='new', action='write', iostat=iostat)
            if (iostat == 0) then
                w%is_open = .true.
            else
                call fail(w, ERR_FILE_OPEN)
            end if
        end if

        ierr = w%ierr
    end subroutine json_open

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Finish the document, close the file and return the writer's first error.
    !|
    !| The document must be complete: a root value written and every container closed,
    !| otherwise the file is deleted and `ERR_JSON_SEQUENCE` returned. A final line feed is
    !| written, the file is closed, and its size on disk is checked against the bytes written.
    !| Whatever fails, the file is gone afterwards, so `ERR_OK` is the only result that leaves
    !| one.
    subroutine json_close(w, ierr)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int32), intent(out) :: ierr
            !! The writer's first error (`w%ierr`), `ERR_OK` when the file is complete

        logical(c_bool) :: may_act
        integer(int32) :: iostat
        integer(int64) :: file_size

        call start_call(w, may_act)
        ! the root is written only when its last container closes, so it also means depth 0
        if (may_act .and. .not. w%root_written) call fail(w, ERR_JSON_SEQUENCE)
        call put(w, new_line('a'))
        call flush_buffer(w)

        if (is_ok(w%ierr)) then
            close (w%unit, iostat=iostat)
            w%is_open = .false.
            file_size = -1
            if (iostat == 0) inquire (file=w%filename, size=file_size, iostat=iostat)
            if (iostat /= 0 .or. file_size /= w%bytes_written) then
                call fail(w, ERR_WRITE_DATA)
                call delete_file_by_name(w%filename)
            end if
        end if

        if (allocated(w%buffer)) deallocate (w%buffer)
        if (allocated(w%container_kinds)) deallocate (w%container_kinds)
        ierr = w%ierr
    end subroutine json_close

    ! ============================================================================================
    ! Containers and null
    ! ============================================================================================

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Open an object: a member of the enclosing object when `key` is given, otherwise an
    !| element of the enclosing array or the root.
    subroutine json_begin_object(w, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in), optional :: key
            !! Key of the object in the enclosing object; required there and refused anywhere else

        call begin_container(w, CONTAINER_OBJECT, key)
    end subroutine json_begin_object

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Open an array: a member of the enclosing object when `key` is given, otherwise an
    !| element of the enclosing array or the root.
    subroutine json_begin_array(w, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in), optional :: key
            !! Key of the array in the enclosing object; required there and refused anywhere else

        call begin_container(w, CONTAINER_ARRAY, key)
    end subroutine json_begin_array

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Close the innermost open container, which must be an object.
    subroutine json_end_object(w)
        type(json_writer), intent(inout) :: w
            !! The writer

        call end_container(w, CONTAINER_OBJECT)
    end subroutine json_end_object

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Close the innermost open container, which must be an array.
    subroutine json_end_array(w)
        type(json_writer), intent(inout) :: w
            !! The writer

        call end_container(w, CONTAINER_ARRAY)
    end subroutine json_end_array

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write `null`: a member of the enclosing object when `key` is given, otherwise an element
    !| of the enclosing array or the root.
    subroutine json_null(w, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in), optional :: key
            !! Key in the enclosing object; required there and refused anywhere else

        logical(c_bool) :: accepted

        call start_value(w, accepted, key)
        if (.not. accepted) return
        call put(w, 'null')
        call finish_value(w)
    end subroutine json_null

    ! ============================================================================================
    ! Checks a caller may run before opening a file
    ! ============================================================================================

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether `text`, without its trailing blanks, is valid UTF-8 -- the test the writer applies
    !| to every key and string value, see [[f42_serde_json_serialize(module)]].
    pure function json_is_valid_utf8(text) result(is_valid)
        character(len=*), intent(in) :: text
            !! The string to test
        logical(c_bool) :: is_valid
            !! `.true.` when every byte from 0x80 up belongs to a well-formed UTF-8 sequence

        integer(int64) :: i_byte, text_length
        integer(int32) :: sequence_length

        is_valid = .false.
        text_length = len_trim(text, kind=int64)
        i_byte = 1
        do while (i_byte <= text_length)
            if (byte_at(text, i_byte) < 128) then
                i_byte = i_byte + 1
            else
                sequence_length = utf8_sequence_length(text, i_byte, text_length)
                if (sequence_length == 0) return
                i_byte = i_byte + sequence_length
            end if
        end do
        is_valid = .true.
    end function json_is_valid_utf8

    ! ============================================================================================
    ! The generics' specifics: a scalar or a rank-1 array of any supported type
    ! ============================================================================================

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write a scalar member; see [[f42_serde_json_serialize(module):json_member(interface)]]
    subroutine member_value(w, key, value)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in) :: key
            !! Key in the enclosing object
        class(*), intent(in) :: value
            !! Value to write, of a supported type
        call write_value(w, value, key)
    end subroutine member_value

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write an array member; see [[f42_serde_json_serialize(module):json_member(interface)]]
    subroutine member_values(w, key, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in) :: key
            !! Key in the enclosing object
        class(*), intent(in) :: values(:)
            !! Values to write, of a supported type
        call write_values(w, values, key)
    end subroutine member_values

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write a scalar element; see [[f42_serde_json_serialize(module):json_element(interface)]]
    subroutine element_value(w, value)
        type(json_writer), intent(inout) :: w
            !! The writer
        class(*), intent(in) :: value
            !! Value to write, of a supported type
        call write_value(w, value)
    end subroutine element_value

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write an array element; see [[f42_serde_json_serialize(module):json_element(interface)]]
    subroutine element_values(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        class(*), intent(in) :: values(:)
            !! Values to write, of a supported type
        call write_values(w, values)
    end subroutine element_values

    ! ============================================================================================
    ! One value, with an optional key: what both generics forward to
    ! ============================================================================================

    !> Write a scalar as the next value, dispatching on its type.
    subroutine write_value(w, value, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        class(*), intent(in) :: value
            !! Value to write
        character(len=*), intent(in), optional :: key
            !! Key in the enclosing object

        logical(c_bool) :: accepted
        character(len=REAL_FIELD_WIDTH) :: field

        call start_value(w, accepted, key)
        if (.not. accepted) return
        select type (value)
        type is (integer(int32))
            write (field, INTEGER_FORMAT) value
            call put(w, trim(field))
        type is (integer(int64))
            write (field, INTEGER_FORMAT) value
            call put(w, trim(field))
        type is (real(real64))
            if (.not. ieee_is_finite(value)) then
                call fail(w, ERR_NAN_INF)
                return
            end if
            write (field, REAL_FORMAT) value
            call put_real_field(w, value, field)
        type is (logical(c_bool))
            call put_boolean(w, value)
        type is (logical)
            call put_boolean(w, logical(value, c_bool))
        type is (character(len=*))
            call put_string(w, value)
        class default
            call fail(w, ERR_INVALID_INPUT)
            return
        end select
        call finish_value(w)
    end subroutine write_value

    !> Write a rank-1 array as the next value, dispatching on its type; an empty one is `[]`.
    !|
    !| Each branch hands the whole array to a routine of its type. Indexing the array inside
    !| `select type` would be shorter, but nvfortran (still 26.9) then misplaces the elements of
    !| an array section, while the whole array passed on arrives intact on every compiler.
    subroutine write_values(w, values, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        class(*), intent(in) :: values(:)
            !! Values to write
        character(len=*), intent(in), optional :: key
            !! Key in the enclosing object

        logical(c_bool) :: accepted

        call start_value(w, accepted, key)
        if (.not. accepted) return
        call put(w, '[')
        ! An empty array is `[]` whatever its type, so the type is only asked for when there
        ! are elements. That is also what keeps nvfortran 26.9 right: it sends a zero-size typed
        ! array constructor such as `[integer ::]` to `class default`.
        if (size(values, kind=int64) > 0) then
            select type (values)
            type is (integer(int32))
                call put_int32_list(w, values)
            type is (integer(int64))
                call put_int64_list(w, values)
            type is (real(real64))
                call put_real_list(w, values)
            type is (logical(c_bool))
                call put_c_bool_list(w, values)
            type is (logical)
                call put_logical_list(w, values)
            type is (character(len=*))
                call put_string_list(w, values)
            class default
                call fail(w, ERR_INVALID_INPUT)
                return
            end select
        end if
        call put(w, ']')
        call finish_value(w)
    end subroutine write_values

    !> Append `int32` values, comma-separated, formatting a block of them per write.
    subroutine put_int32_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int32), intent(in) :: values(:)
            !! The values

        character(len=INTEGER_FIELD_WIDTH*NUMBER_BLOCK) :: block_text
        integer(int64) :: block_start, block_end, i_value

        do block_start = 1_int64, size(values, kind=int64), int(NUMBER_BLOCK, int64)
            block_end = min(block_start + NUMBER_BLOCK - 1, size(values, kind=int64))
            write (block_text, INTEGER_FORMAT) (values(i_value), i_value = block_start, block_end)
            if (block_start > 1) call put(w, ',')
            call put(w, trim(block_text))
        end do
    end subroutine put_int32_list

    !> Append `int64` values, comma-separated, formatting a block of them per write.
    subroutine put_int64_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int64), intent(in) :: values(:)
            !! The values

        character(len=INTEGER_FIELD_WIDTH*NUMBER_BLOCK) :: block_text
        integer(int64) :: block_start, block_end, i_value

        do block_start = 1_int64, size(values, kind=int64), int(NUMBER_BLOCK, int64)
            block_end = min(block_start + NUMBER_BLOCK - 1, size(values, kind=int64))
            write (block_text, INTEGER_FORMAT) (values(i_value), i_value = block_start, block_end)
            if (block_start > 1) call put(w, ',')
            call put(w, trim(block_text))
        end do
    end subroutine put_int64_list

    !> Append finite reals, comma-separated, formatting a block of them per write.
    subroutine put_real_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        real(real64), intent(in) :: values(:)
            !! The values; all must be finite

        character(len=REAL_FIELD_WIDTH*NUMBER_BLOCK) :: block_text
        integer(int64) :: block_start, block_end, i_value, field_start

        do block_start = 1_int64, size(values, kind=int64), int(NUMBER_BLOCK, int64)
            block_end = min(block_start + NUMBER_BLOCK - 1, size(values, kind=int64))
            if (.not. all(ieee_is_finite(values(block_start:block_end)))) then
                call fail(w, ERR_NAN_INF)
                return
            end if
            ! an implied do, not a section: a section of a strided array is copied first by ifx
            write (block_text, REAL_FORMAT) (values(i_value), i_value = block_start, block_end)
            do i_value = block_start, block_end
                if (i_value > 1) call put(w, ',')
                field_start = (i_value - block_start)*REAL_FIELD_WIDTH
                call put_real_field(w, values(i_value), block_text(field_start + 1:field_start + REAL_FIELD_WIDTH))
            end do
        end do
    end subroutine put_real_list

    !> Append `logical(c_bool)` values as `true`/`false`, comma-separated.
    subroutine put_c_bool_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        logical(c_bool), intent(in) :: values(:)
            !! The values

        integer(int64) :: i_value

        do i_value = 1_int64, size(values, kind=int64)
            if (i_value > 1) call put(w, ',')
            call put_boolean(w, values(i_value))
        end do
    end subroutine put_c_bool_list

    !> Append default `logical` values as `true`/`false`, comma-separated.
    subroutine put_logical_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        logical, intent(in) :: values(:)
            !! The values

        integer(int64) :: i_value

        do i_value = 1_int64, size(values, kind=int64)
            if (i_value > 1) call put(w, ',')
            call put_boolean(w, logical(values(i_value), c_bool))
        end do
    end subroutine put_logical_list

    !> Append strings, comma-separated.
    subroutine put_string_list(w, values)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in) :: values(:)
            !! The strings

        integer(int64) :: i_value

        do i_value = 1_int64, size(values, kind=int64)
            if (i_value > 1) call put(w, ',')
            call put_string(w, values(i_value))
        end do
    end subroutine put_string_list


    ! ============================================================================================
    ! The state machine
    ! ============================================================================================

    !> Count a call, and tell whether it may act: the writer is open and has not failed. A call
    !| on a writer that is not open fails it with `ERR_JSON_SEQUENCE`.
    subroutine start_call(w, may_act)
        type(json_writer), intent(inout) :: w
            !! The writer
        logical(c_bool), intent(out) :: may_act
            !! Whether the call may go on

        w%call_count = w%call_count + 1
        may_act = .false.
        ! a failed writer is never open: `fail` closes it
        if (.not. w%is_open) then
            call fail(w, ERR_JSON_SEQUENCE)
            return
        end if
        may_act = .true.
    end subroutine start_call

    !> Count the call and check that a value may start here, with or without a key; if it may,
    !| write the comma before it and its key.
    subroutine start_value(w, accepted, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        logical(c_bool), intent(out) :: accepted
            !! Whether the value may be written now
        character(len=*), intent(in), optional :: key
            !! Key in the enclosing object

        logical(c_bool) :: is_misplaced

        call start_call(w, accepted)
        if (.not. accepted) return

        if (w%depth == 0) then
            is_misplaced = present(key) .or. w%root_written
        else
            is_misplaced = (w%container_kinds(w%depth) == CONTAINER_OBJECT) .neqv. present(key)
        end if
        if (is_misplaced) then
            call fail(w, ERR_JSON_SEQUENCE)
            accepted = .false.
            return
        end if

        if (.not. w%container_is_empty) call put(w, ',')
        if (present(key)) then
            call put_string(w, key)
            call put(w, ':')
        end if
        w%container_is_empty = .false.
        accepted = is_ok(w%ierr)
    end subroutine start_value

    !> Record that a value is complete; at the top level it is the root.
    subroutine finish_value(w)
        type(json_writer), intent(inout) :: w
            !! The writer

        if (w%depth == 0 .and. is_ok(w%ierr)) w%root_written = .true.
    end subroutine finish_value

    !> Start a container as the next value and push it on the stack, growing the stack when full.
    subroutine begin_container(w, container_kind, key)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int8), intent(in) :: container_kind
            !! `CONTAINER_OBJECT` or `CONTAINER_ARRAY`
        character(len=*), intent(in), optional :: key
            !! Key in the enclosing object

        logical(c_bool) :: accepted
        integer(int8), allocatable :: grown_kinds(:)
        integer(int32) :: alloc_stat

        call start_value(w, accepted, key)
        if (.not. accepted) return

        if (w%depth == size(w%container_kinds, kind=int32)) then
            allocate (grown_kinds(2*w%depth), stat=alloc_stat)
            if (alloc_stat /= 0) then
                call fail(w, ERR_ALLOC_FAIL)
                return
            end if
            grown_kinds(1:w%depth) = w%container_kinds(1:w%depth)
            call move_alloc(grown_kinds, w%container_kinds)
        end if
        w%depth = w%depth + 1
        w%container_kinds(w%depth) = container_kind
        w%container_is_empty = .true.
        call put(w, OPENING_BRACKETS(container_kind:container_kind))
    end subroutine begin_container

    !> Count the call, check that the innermost open container is of the kind given, and close it.
    subroutine end_container(w, container_kind)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int8), intent(in) :: container_kind
            !! `CONTAINER_OBJECT` or `CONTAINER_ARRAY`

        logical(c_bool) :: may_act

        call start_call(w, may_act)
        if (.not. may_act) return
        if (w%depth == 0) then
            call fail(w, ERR_JSON_SEQUENCE)
            return
        end if
        if (w%container_kinds(w%depth) /= container_kind) then
            call fail(w, ERR_JSON_SEQUENCE)
            return
        end if

        call put(w, CLOSING_BRACKETS(container_kind:container_kind))
        w%depth = w%depth - 1
        w%container_is_empty = .false.
        call finish_value(w)
    end subroutine end_container

    !> Record the first error and the call that caused it, and delete the file at once.
    subroutine fail(w, error_code)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int32), intent(in) :: error_code
            !! The `ERR_*` code

        integer(int32) :: iostat

        if (is_err(w%ierr)) return
        call set_err(w%ierr, error_code)
        w%failing_call = w%call_count
        if (w%is_open) then
            close (w%unit, status='delete', iostat=iostat)
            w%is_open = .false.
            ! A failed close still disconnects the unit, so the file is reached by its name.
            if (iostat /= 0) call delete_file_by_name(w%filename)
        end if
    end subroutine fail

    !> Delete a file by its name, if it exists; the last resort after a failed close.
    subroutine delete_file_by_name(filename)
        character(len=*), intent(in) :: filename
            !! Path of the file

        integer(int32) :: unit, iostat
        ! default kind, the only one the standard allows for INQUIRE's EXIST=
        logical :: file_exists

        inquire (file=filename, exist=file_exists, iostat=iostat)
        if (iostat /= 0 .or. .not. file_exists) return
        open (newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat == 0) close (unit, status='delete', iostat=iostat)
    end subroutine delete_file_by_name

    ! ============================================================================================
    ! Output
    ! ============================================================================================

    !> Append text to the buffer, writing the buffer out first when the text does not fit.
    subroutine put(w, text)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in) :: text
            !! Text to append

        integer(int64) :: text_length
        integer(int32) :: iostat

        if (is_err(w%ierr)) return
        text_length = len(text, kind=int64)
        if (w%buffer_used + text_length > BUFFER_CAPACITY) then
            call flush_buffer(w)
            if (is_err(w%ierr)) return
            if (text_length > BUFFER_CAPACITY) then
                write (w%unit, iostat=iostat) text
                if (iostat /= 0) then
                    call fail(w, ERR_WRITE_DATA)
                    return
                end if
                w%bytes_written = w%bytes_written + text_length
                return
            end if
        end if
        w%buffer(w%buffer_used + 1:w%buffer_used + text_length) = text
        w%buffer_used = w%buffer_used + int(text_length, int32)   ! at most BUFFER_CAPACITY here
    end subroutine put

    !> Write the buffer to the file and empty it.
    subroutine flush_buffer(w)
        type(json_writer), intent(inout) :: w
            !! The writer

        integer(int32) :: iostat

        if (is_err(w%ierr) .or. w%buffer_used == 0) return
        write (w%unit, iostat=iostat) w%buffer(1:w%buffer_used)
        if (iostat /= 0) then
            call fail(w, ERR_WRITE_DATA)
            return
        end if
        w%bytes_written = w%bytes_written + w%buffer_used
        w%buffer_used = 0
    end subroutine flush_buffer

    !> Append `true` or `false`.
    subroutine put_boolean(w, value)
        type(json_writer), intent(inout) :: w
            !! The writer
        logical(c_bool), intent(in) :: value
            !! The value

        if (value) then
            call put(w, 'true')
        else
            call put(w, 'false')
        end if
    end subroutine put_boolean

    !> Append one finite real, given its `ES24.16E3` field; a zero is recognised by its IEEE class.
    !|
    !| Not by comparing with zero: where denormals are flushed (ifx from `-O1` on) a subnormal
    !| compares equal to zero, and `sign` drops the sign of `-0.0` on ifx without
    !| `-assume minus0`. The class is right on every supported compiler and profile.
    subroutine put_real_field(w, value, field)
        type(json_writer), intent(inout) :: w
            !! The writer
        real(real64), intent(in) :: value
            !! The value, finite
        character(len=REAL_FIELD_WIDTH), intent(in) :: field
            !! `value` formatted as `ES24.16E3`

        type(ieee_class_type) :: value_class

        value_class = ieee_class(value)
        if (value_class == ieee_negative_zero) then
            call put(w, '-0.0')
        else if (value_class == ieee_positive_zero) then
            call put(w, '0.0')
        else
            call put(w, field(verify(field, ' '):))
        end if
    end subroutine put_real_field

    !> Append a string in quotes: trailing blanks removed, escapes applied, UTF-8 validated.
    !|
    !| Bytes that need no escape are appended in runs, one `put` per run.
    subroutine put_string(w, text)
        type(json_writer), intent(inout) :: w
            !! The writer
        character(len=*), intent(in) :: text
            !! The string

        integer(int64) :: text_length, i_byte, run_start
        integer(int32) :: byte_value, sequence_length

        call put(w, '"')
        text_length = len_trim(text, kind=int64)
        run_start = 1
        i_byte = 1
        do while (i_byte <= text_length)
            byte_value = byte_at(text, i_byte)
            if (byte_value >= 128) then
                sequence_length = utf8_sequence_length(text, i_byte, text_length)
                if (sequence_length == 0) then
                    call fail(w, ERR_INVALID_UTF8)
                    return
                end if
                i_byte = i_byte + sequence_length
            else if (byte_value >= 32 .and. byte_value /= 34 .and. byte_value /= 92) then
                i_byte = i_byte + 1
            else
                if (i_byte > run_start) call put(w, text(run_start:i_byte - 1))
                call put_escape(w, byte_value)
                i_byte = i_byte + 1
                run_start = i_byte
            end if
        end do
        if (text_length >= run_start) call put(w, text(run_start:text_length))
        call put(w, '"')
    end subroutine put_string

    !> Append the escape sequence of `"`, `\` or a control character.
    subroutine put_escape(w, byte_value)
        type(json_writer), intent(inout) :: w
            !! The writer
        integer(int32), intent(in) :: byte_value
            !! The byte, 34, 92 or below 32

        select case (byte_value)
        case (34)
            call put(w, BACKSLASH//'"')
        case (92)
            call put(w, BACKSLASH//BACKSLASH)
        case (8)
            call put(w, BACKSLASH//'b')
        case (9)
            call put(w, BACKSLASH//'t')
        case (10)
            call put(w, BACKSLASH//'n')
        case (12)
            call put(w, BACKSLASH//'f')
        case (13)
            call put(w, BACKSLASH//'r')
        case default
            call put(w, BACKSLASH//'u00'//HEX_DIGITS(byte_value/16 + 1:byte_value/16 + 1) &
                     //HEX_DIGITS(mod(byte_value, 16) + 1:mod(byte_value, 16) + 1))
        end select
    end subroutine put_escape

    ! ============================================================================================
    ! UTF-8
    ! ============================================================================================

    !> The byte at a position of a string, as a value from 0 to 255.
    pure function byte_at(text, position) result(byte_value)
        character(len=*), intent(in) :: text
            !! The string
        integer(int64), intent(in) :: position
            !! Position of the byte
        integer(int32) :: byte_value
            !! Its value

        byte_value = iand(ichar(text(position:position), kind=int32), 255_int32)
    end function byte_at

    !> Length of the well-formed UTF-8 sequence starting at a byte from 0x80 up, 0 if it is not one.
    !|
    !| Follows the table of well-formed byte sequences of the Unicode standard (section 3.9): the
    !| lead byte fixes the length and narrows the range of the second byte, which excludes the
    !| overlong forms, the surrogates and everything above U+10FFFF.
    pure function utf8_sequence_length(text, start, text_end) result(sequence_length)
        character(len=*), intent(in) :: text
            !! The string
        integer(int64), intent(in) :: start
            !! Position of the lead byte
        integer(int64), intent(in) :: text_end
            !! Position of the last byte that belongs to the string
        integer(int32) :: sequence_length
            !! Bytes in the sequence (2 to 4), or 0 when it is malformed

        integer(int32) :: lead_byte, expected_length, second_min, second_max, continuation
        integer(int64) :: i_continuation

        sequence_length = 0
        lead_byte = byte_at(text, start)
        second_min = 128
        second_max = 191
        select case (lead_byte)
        case (194:223)
            expected_length = 2
        case (224:239)
            expected_length = 3
            if (lead_byte == 224) second_min = 160
            if (lead_byte == 237) second_max = 159
        case (240:244)
            expected_length = 4
            if (lead_byte == 240) second_min = 144
            if (lead_byte == 244) second_max = 143
        case default
            return
        end select
        if (start + expected_length - 1 > text_end) return

        continuation = byte_at(text, start + 1)
        if (continuation < second_min .or. continuation > second_max) return
        do i_continuation = start + 2, start + expected_length - 1
            continuation = byte_at(text, i_continuation)
            if (continuation < 128 .or. continuation > 191) return
        end do
        sequence_length = expected_length
    end function utf8_sequence_length

end module f42_serde_json_serialize
