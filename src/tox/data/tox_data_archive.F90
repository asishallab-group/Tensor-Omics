#include <src/macros.h>

!> Zip-archive backed persistence for TensorOmics data sets.
!|
!| Wraps libzip (via C bindings declared in the interface block below) to create/extract zip
!| archives whose members are [[tox_data_read_write(module)]]-serialized arrays, indexed by a
!| plain-text `manifest.txt` mapping logical keys (e.g. `gene_ids`, `expression`) to member
!| filenames. [[tox_data_archive(module):save_tox_data(subroutine)]] and
!| [[tox_data_archive(module):read_tox_data(subroutine)]] are the standard entry points for the
!| fixed TensorOmics data set schema; `create_zip_archive`/`extract_zip_archive` and the
!| `*_manifest*` routines below are the generic key/filename building blocks they are built on.
module tox_data_archive
    use safeguard
    use iso_c_binding, only: c_ptr, c_char, c_int, c_int64_t, c_size_t, c_signed_char, c_f_pointer, c_loc, c_associated, c_null_char, c_null_ptr
    use tox_data_read_write
    use tox_errors, only: set_ok, set_err_once, is_err, ERR_FILE_OPEN, ERR_ALLOC_FAIL, ERR_FILE_ADD, set_err
    use tox_errors, only: ERR_FILE_CLOSE, ERR_FILE_EXTRACT, ERR_INVALID_INPUT
    use tox_errors, only: ERR_POINTER_NULL, ERR_WRITE_DATA, ERR_READ_DATA
    use iso_fortran_env, only: real64, int32, iostat_end
    use config, only: DEBUG
    implicit none

    ! libzip constants
    integer(c_int), parameter :: ZIP_CREATE = 1
    integer(c_int), parameter :: ZIP_EXCLUSIVE = 3
    integer(c_int), parameter :: ZIP_READ_ONLY = 0
    integer(c_int), parameter :: ZIP_FILE_OVERWRITE = 8192
    integer(c_int), parameter :: ZIP_COMPRESSION_STORE = 0
    integer(c_int), parameter :: ZIP_COMPRESSION_DEFLATE = 8

    ! Constants for data types
    integer(int32), parameter :: DATA_TYPE_FILE = 1
    integer(int32), parameter :: DATA_TYPE_STRING = 2

    ! libzip interface definitions
    interface
        function malloc(size) bind(C, name="malloc")
            use iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr) :: malloc
            integer(c_size_t), value :: size
        end function malloc

        subroutine free(ptr) bind(C, name="free")
            use iso_c_binding, only: c_ptr
            type(c_ptr), value :: ptr
        end subroutine free

        function zip_open(path, flags, errorp) bind(C, name="zip_open")
            use iso_c_binding, only: c_ptr, c_char, c_int
            type(c_ptr) :: zip_open
            character(kind=c_char), dimension(*) :: path
            integer(c_int), value :: flags
            integer(c_int), intent(out) :: errorp
        end function zip_open

        function zip_close(archive) bind(C, name="zip_close")
            use iso_c_binding, only: c_ptr, c_int
            integer(c_int) :: zip_close
            type(c_ptr), value :: archive
        end function zip_close

        function zip_file_add(archive, name, source, flags) bind(C, name="zip_file_add")
            use iso_c_binding, only: c_int64_t, c_ptr, c_char, c_int
            integer(c_int64_t) :: zip_file_add
            type(c_ptr), value :: archive
            character(kind=c_char), dimension(*) :: name
            type(c_ptr), value :: source
            integer(c_int), value :: flags
        end function zip_file_add

        function zip_source_buffer(archive, data, len, freep) bind(C, name="zip_source_buffer")
            use iso_c_binding, only: c_ptr, c_size_t, c_int
            type(c_ptr) :: zip_source_buffer
            type(c_ptr), value :: archive
            type(c_ptr), value :: data
            integer(c_size_t), value :: len
            integer(c_int), value :: freep
        end function zip_source_buffer

        function zip_set_file_compression(archive, index, comp, flags) bind(C, name="zip_set_file_compression")
            use iso_c_binding, only: c_int, c_ptr, c_int64_t
            integer(c_int) :: zip_set_file_compression
            type(c_ptr), value :: archive
            integer(c_int64_t), value :: index
            integer(c_int), value :: comp
            integer(c_int), value :: flags
        end function zip_set_file_compression

        function zip_fopen(archive, fname, flags) bind(C, name="zip_fopen")
            use iso_c_binding, only: c_ptr, c_char, c_int
            type(c_ptr) :: zip_fopen
            type(c_ptr), value :: archive
            character(kind=c_char), dimension(*) :: fname
            integer(c_int), value :: flags
        end function zip_fopen

        function zip_fread(file, buf, nbytes) bind(C, name="zip_fread")
            use iso_c_binding, only: c_int64_t, c_ptr, c_size_t
            integer(c_int64_t) :: zip_fread
            type(c_ptr), value :: file
            type(c_ptr), value :: buf
            integer(c_size_t), value :: nbytes
        end function zip_fread

        function zip_fclose(file) bind(C, name="zip_fclose")
            use iso_c_binding, only: c_int, c_ptr
            integer(c_int) :: zip_fclose
            type(c_ptr), value :: file
        end function zip_fclose

        function zip_get_num_entries(archive, flags) bind(C, name="zip_get_num_entries")
            use iso_c_binding, only: c_int64_t, c_ptr, c_int
            integer(c_int64_t) :: zip_get_num_entries
            type(c_ptr), value :: archive
            integer(c_int), value :: flags
        end function zip_get_num_entries

        function zip_get_name(archive, index, flags) bind(C, name="zip_get_name")
            use iso_c_binding, only: c_ptr, c_int64_t, c_int
            type(c_ptr) :: zip_get_name
            type(c_ptr), value :: archive
            integer(c_int64_t), value :: index
            integer(c_int), value :: flags
        end function zip_get_name

        subroutine zip_source_free(source) bind(C, name="zip_source_free")
            use iso_c_binding, only: c_ptr
            type(c_ptr), value :: source
        end subroutine zip_source_free
    end interface

contains

    !> AUTHOR_AARON_SCHROEDER
    !| Creates a zip archive with generic file lists
    subroutine create_zip_archive(zip_filename, keys, filenames, ierr)
        character(len=*), intent(in) :: zip_filename
            !! Name of the zip file to create
        character(len=*), intent(in) :: keys(:)
            !! Array of keys for manifest entries
        character(len=*), intent(in) :: filenames(:)
            !! Array of filenames to add to zip
        integer(int32), intent(out) :: ierr
            !! Error code

        type(c_ptr) :: zip_handle
        integer(c_int) :: error
        character(len=:), allocatable :: manifest_filename
        integer(int32) :: i

        call set_ok(ierr)
        call set_ok(error)

        ! Validate input arrays
        if (size(keys) /= size(filenames)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if (DEBUG) print *, "Error: keys and filenames arrays must have same size"
            return
        end if

        if (len_trim(zip_filename) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if (DEBUG) print *, "Error: zip_filename cannot be empty"
            return
        end if

        ! ZIP_EXCLUSIVE makes zip_open fail (rather than truncate/overwrite) if zip_filename
        ! already exists, so a stale archive is never silently clobbered by a repeated run.
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_EXCLUSIVE, error)
        if (is_err(error)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (error == 10) print *, "Error opening ZIP file for writing: File already exists"
            return
        end if

        ! Add all data files
        do i = 1, size(filenames)
            if (len_trim(filenames(i)) > 0) then
                call add_data_to_zip(zip_handle, filenames(i), filenames(i), DATA_TYPE_FILE, ierr)
                if (is_err(ierr)) then
                    !CLAUDE: On error the partially-written zip_filename is never deleted from disk. Since it was opened with ZIP_EXCLUSIVE, a retry of create_zip_archive with the same name will now fail with "file already exists".
                    error = zip_close(zip_handle)
                    return
                end if
            end if
        end do

        ! Add manifest
        manifest_filename = "manifest.txt"
        call write_manifest(keys, filenames, manifest_filename, ierr)
        if (is_err(ierr)) then
            error = zip_close(zip_handle)
            return
        end if

        call add_data_to_zip(zip_handle, manifest_filename, manifest_filename, DATA_TYPE_FILE, ierr)
        if (is_err(ierr)) then
            error = zip_close(zip_handle)
            return
        end if

        call delete_file(manifest_filename, ierr)
        if (is_err(ierr)) then
            error = zip_close(zip_handle)
            return
        end if

        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (is_err(error)) then
            call set_err_once(ierr, error)
            if (DEBUG) print *, "Error closing ZIP file: ", error
        else if (is_ok(ierr)) then
            if (DEBUG) print *, "ZIP archive created successfully: ", trim(zip_filename)
        end if
    end subroutine create_zip_archive

    !> AUTHOR_AARON_SCHROEDER
    !| Extract a zip archive and return all key-value pairs from manifest
    subroutine extract_zip_archive(zip_filename, keys, filenames, ierr)
        use tox_conversions, only: int32_as_c_int64
        character(len=*), intent(in) :: zip_filename
            !! Zip file to read
        character(len=:), allocatable, intent(out) :: keys(:)
            !! Array of keys from manifest
        character(len=:), allocatable, intent(out) :: filenames(:)
            !! Array of filenames from manifest
        integer(int32), intent(out) :: ierr
            !! Error code

        type(c_ptr) :: zip_handle
        integer(c_int) :: error
        integer(c_int64_t) :: i, num_entries
        character(len=:), allocatable :: filename
        logical :: file_exists
        integer(int32) :: i_fortran
        integer(c_int64_t) :: i_c

        ! Initialize outputs
        M_ALLOCATE(character(len=256) :: keys(0))
        M_ALLOCATE(character(len=256) :: filenames(0))

        call set_ok(ierr)
        call set_ok(error)

        ! Check if file exists
        inquire (file=zip_filename, exist=file_exists)
        if (.not. file_exists) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (DEBUG) print *, "ZIP file does not exist: ", trim(zip_filename)
            return
        end if

        ! Open ZIP archive
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_READ_ONLY, error)
        if (error /= 0 .or. .not. c_associated(zip_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (DEBUG) print *, "Error opening ZIP file for reading: ", error
            return
        end if

        ! Extract all files
        num_entries = zip_get_num_entries(zip_handle, 0)
        do i_fortran = 0, int(num_entries - 1, int32)
            call int32_as_c_int64(i_fortran, i_c)
            call get_zip_entry_name(zip_handle, i_c, filename, ierr)
            if (is_err(ierr)) then
                error = zip_close(zip_handle)
                return
            end if

            ! Reject entry names that could escape the extraction directory (zip-slip)
            if (index(filename, "..") > 0 .or. &
                (len_trim(filename) > 0 .and. (filename(1:1) == "/" .or. filename(1:1) == "\"))) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if (DEBUG) print *, "Error: unsafe ZIP entry name rejected: ", trim(filename)
                error = zip_close(zip_handle)
                return
            end if

            call extract_file_from_zip(zip_handle, filename, ierr)
            if (is_err(ierr)) then
                error = zip_close(zip_handle)
                return
            end if
        end do

        ! Extract and parse manifest using generic version
        call extract_and_parse_manifest(zip_handle, keys, filenames, ierr)
        if (is_err(ierr)) then
            error = zip_close(zip_handle)
            return
        end if

        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (is_err(error)) then
            if (DEBUG) print *, "Error closing ZIP file: ", error
            call set_err_once(ierr, ERR_FILE_CLOSE)
            return
        end if

        if (DEBUG) print *, "ZIP archive extracted successfully: ", trim(zip_filename)
    end subroutine extract_zip_archive

    !> AUTHOR_AARON_SCHROEDER
    !| Delete a file from the disk
    subroutine delete_file(filename, ierr)
        character(len=*), intent(in) :: filename
            !! File to delete
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: unit, iostat

        call set_ok(iostat)
        call set_ok(ierr)

        open (newunit=unit, file=filename, iostat=iostat, status='old')
        if (is_ok(iostat)) then
            close (unit, status='delete')
        else
            call set_err_once(ierr, ERR_FILE_OPEN)
        end if
    end subroutine delete_file

    !> AUTHOR_AARON_SCHROEDER
    !| Helper function to get the name of a ZIP entry
    subroutine get_zip_entry_name(zip_handle, entry_index, entry_name, ierr)
        use tox_conversions, only: c_char_1d_as_string
        type(c_ptr), intent(in) :: zip_handle
            !! Zip file connection
        integer(c_int64_t), intent(in) :: entry_index
            !! Index of the entry
        character(len=:), allocatable, intent(out) :: entry_name
            !! Name of the entry
        integer(int32), intent(out) :: ierr
            !! Error code

        character(kind=c_char, len=1), pointer :: name_ptr(:)
        integer(int32), parameter :: MAX_NAME_LENGTH = 4096  ! Reasonable maximum
        ! zip_get_name returns a pointer to a NUL-terminated C string of unknown length; we don't
        ! have a real bound from libzip, so we map a fixed-size window over it and rely on
        ! c_char_1d_as_string to stop at the first NUL byte. MAX_NAME_LENGTH must stay >= the
        ! longest entry name any archive we read can contain.
        call c_f_pointer(zip_get_name(zip_handle, entry_index, 0), name_ptr, [MAX_NAME_LENGTH])
        call c_char_1d_as_string(name_ptr, entry_name, ierr)
    end subroutine get_zip_entry_name

    !> AUTHOR_AARON_SCHROEDER
    !| Extracts a single named entry from an open ZIP archive to a file of the same name on disk,
    !| streaming it through a fixed-size buffer rather than reading it into memory whole.
    subroutine extract_file_from_zip(zip_handle, filename, ierr)
        use tox_conversions, only: int32_as_c_size
        type(c_ptr), intent(in) :: zip_handle
            !! Zip connection
        character(len=*), intent(in) :: filename
            !! Name of the file to extract
        integer(int32), intent(out) :: ierr
            !! Error code

        type(c_ptr) :: file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: bytes_read
        integer(int32) :: unit, iostat
        integer(int32), parameter :: CHUNK_SIZE = 4096
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        integer(c_size_t) :: chunk_size_c

        call set_ok(ierr)

        ! Open file in ZIP
        file_handle = zip_fopen(zip_handle, trim(filename)//c_null_char, 0)
        if (.not. c_associated(file_handle)) then
            call set_err_once(ierr, ERR_FILE_EXTRACT)
            if (DEBUG) print *, "Error opening file in ZIP: ", trim(filename)
            return
        end if

        ! Open output file
        open (newunit=unit, file=trim(filename), access='stream', form='unformatted', &
              iostat=iostat, status='replace', action='write')
        if (iostat /= 0) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (DEBUG) print *, "Error creating file: ", trim(filename)
            error = zip_fclose(file_handle)
            return
        end if

        ! Read and write in chunks
        M_ALLOCATE(buffer(CHUNK_SIZE))

        call int32_as_c_size(CHUNK_SIZE, chunk_size_c)
        do
            bytes_read = zip_fread(file_handle, c_loc(buffer), chunk_size_c)
            if (bytes_read <= 0) exit
            if (bytes_read > 0) then
                write (unit, iostat=iostat) buffer(1:bytes_read)
                if (iostat /= 0) then
                    call set_err_once(ierr, ERR_WRITE_DATA)
                    if (DEBUG) print *, "Error writing file: ", trim(filename)
                    exit
                end if
            end if
        end do

        ! Clean up
        if (allocated(buffer)) deallocate (buffer)
        close (unit)
        error = zip_fclose(file_handle)

        if (is_err(error)) then
            call set_err_once(ierr, ERR_FILE_CLOSE)
            if (DEBUG) print *, "Error closing file in ZIP: ", trim(filename)
        end if

        if (is_ok(ierr)) then
            if (DEBUG) print *, "Extracted: ", trim(filename)
        end if
    end subroutine extract_file_from_zip

    !> AUTHOR_AARON_SCHROEDER
    !| Adds one entry to an open ZIP archive, sourced either from a file on disk
    !| (`data_type=DATA_TYPE_FILE`, `data_source` is a path) or from an in-memory string
    !| (`data_type=DATA_TYPE_STRING`, `data_source` is the literal content).
    subroutine add_data_to_zip(zip_handle, filename, data_source, data_type, ierr)
        use tox_conversions, only: int32_as_c_size
        implicit none

        ! Arguments
        type(c_ptr), intent(in)    :: zip_handle
            !! Zip connection
        character(len=*), intent(in) :: filename
            !! Filename to add
        character(len=*), intent(in) :: data_source
            !! File path or string content
        integer(int32), intent(in) :: data_type
            !! Type of input
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Locals
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source, c_data
        integer(c_signed_char), pointer :: file_data(:)
        character(kind=c_char), pointer :: string_data(:)
        integer(int32) :: unit, iostat, file_size, i
        integer(c_size_t) :: data_len

        call set_ok(ierr)
        c_data = c_null_ptr

        if (len_trim(data_source) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! Each branch below allocates its payload with C malloc (not Fortran allocate) and then
        ! hands it to zip_source_buffer with freep=1, which tells libzip it now owns the buffer
        ! and must free() it once the source is consumed/closed. Do not deallocate c_data
        ! ourselves after a successful zip_source_buffer call, or it will be double-freed; on the
        ! error paths before that handoff happens, this subroutine is responsible for free()-ing
        ! it manually instead.
        select case (data_type)

        case (DATA_TYPE_FILE)
            ! Open file
            open (newunit=unit, file=data_source, access='stream', form='unformatted', &
                  iostat=iostat, status='old')
            if (iostat /= 0) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                if (DEBUG) print *, "Error opening file: ", trim(data_source)
                return
            end if

            inquire (unit, size=file_size, iostat=iostat)
            if (iostat /= 0 .or. file_size < 0) then
                call set_err_once(ierr, ERR_READ_DATA)
                if (DEBUG) print *, "Error: could not determine size of file: ", trim(data_source)
                close (unit)
                return
            end if
            call int32_as_c_size(file_size, data_len)

            if (file_size == 0) then
                close (unit)
                ! No buffer to hand over, so freep=0: there is nothing for libzip to free.
                source = zip_source_buffer(zip_handle, c_null_ptr, 0_c_size_t, 0)
            else
                c_data = malloc(data_len)
                if (.not. c_associated(c_data)) then
                    call set_err_once(ierr, ERR_POINTER_NULL)
                    close (unit)
                    return
                end if

                call c_f_pointer(c_data, file_data, [file_size])
                read (unit, iostat=iostat) file_data
                close (unit)

                if (iostat /= 0) then
                    call set_err_once(ierr, ERR_READ_DATA)
                    call free(c_data)
                    return
                end if

                source = zip_source_buffer(zip_handle, c_data, data_len, 1)
            end if

        case (DATA_TYPE_STRING)
            data_len = len(data_source)
            c_data = malloc(data_len)
            if (.not. c_associated(c_data)) then
                call set_err_once(ierr, ERR_POINTER_NULL)
                return
            end if

            call c_f_pointer(c_data, string_data, [data_len])
            do i = 1, data_len
                string_data(i) = data_source(i:i)
            end do

            source = zip_source_buffer(zip_handle, c_data, data_len, 1)

        case default
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end select

        if (.not. c_associated(source)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            if (c_associated(c_data)) call free(c_data)
            return
        end if

        ! Add file to ZIP
        index = zip_file_add(zip_handle, trim(filename)//c_null_char, source, ZIP_FILE_OVERWRITE)
        if (index < 0) then
            call set_err_once(ierr, ERR_FILE_ADD)
            call zip_source_free(source)
            return
        end if

        ! Set compression to "store" (no compression)
        error = zip_set_file_compression(zip_handle, index, ZIP_COMPRESSION_STORE, 0)
        if (is_err(error)) then
            if (DEBUG) print *, "Warning: Error setting compression for: ", trim(filename)
        end if

        if (DEBUG) print *, "Added to ZIP: ", trim(filename)

    end subroutine add_data_to_zip

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Write manifest from given key-value pairs
    subroutine write_manifest(keys, filenames, manifest_filename, ierr)
        character(len=*), intent(in) :: keys(:)
            !! Array of keys for manifest entries
        character(len=*), intent(in) :: filenames(:)
            !! Array of filenames for manifest entries
        character(len=*), intent(in) :: manifest_filename
            !! Name of the manifest file (should be manifest.txt)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit, iostat, i_key

        call set_ok(ierr)

        ! Validate input arrays have same size
        if (size(keys) /= size(filenames)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if (DEBUG) print *, "Error: keys and filenames arrays must have same size"
            return
        end if

        ! Open the manifest file for writing
        open (newunit=unit, file=manifest_filename, status='replace', iostat=iostat)
        if (iostat /= 0) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (DEBUG) print *, "Error creating manifest file: ", trim(manifest_filename)
            return
        end if

        write(unit, '(I0)') size(keys)
        write(unit, '(I0)') len(keys)
        write(unit, '(I0)') len(filenames)
        ! Write each key-value pair to the manifest
        do i_key = 1, size(keys)
            if (len_trim(keys(i_key)) > 0 .and. len_trim(filenames(i_key)) > 0) then
                write(unit, '(a)') trim(keys(i_key))
                write(unit, '(a)') trim(filenames(i_key))
            else
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if (DEBUG) print *, "Error: keys and filenames must not be empty strings"
                close (unit)
                return
            end if
        end do

        ! Close the manifest file
        close (unit)

        if (DEBUG) print *, "Manifest created successfully with ", size(keys), " entries"
    end subroutine write_manifest

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Read manifest file and return key-value pairs
    subroutine read_manifest_generic(manifest_filename, keys, values, ierr)
        character(len=*), intent(in) :: manifest_filename
            !! Filename of the manifest (should be manifest.txt)
        character(len=:), allocatable, intent(out) :: keys(:)
            !! Array of keys from manifest
        character(len=:), allocatable, intent(out) :: values(:)
            !! Array of values from manifest
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit, iostat, i_pair, strlen, n_pairs

        call set_ok(ierr)

        ! Open the manifest file for reading
        open (newunit=unit, file=manifest_filename, status='old', iostat=iostat, action='read')
        if (iostat /= 0) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (DEBUG) print *, "Error opening manifest file: ", trim(manifest_filename)
            return
        end if

#define CM_CORRUPTED then; call set_err_once(ierr, ERR_READ_DATA); if (DEBUG) print *, "Error reading manifest file: ", trim(manifest_filename); close(unit); return; end if
#define CM_CHECK_CORRUPTED if (iostat /= 0) CM_CORRUPTED
#define CM_READ(TARGET) read(unit, *, iostat=iostat) TARGET; CM_CHECK_CORRUPTED
#define CM_READ_DIM(TARGET) CM_READ(TARGET); if (TARGET <= 0) CM_CORRUPTED

        CM_READ_DIM(n_pairs)
        CM_READ_DIM(strlen)
        M_ALLOCATE(character(len=strlen) :: keys(n_pairs))
        CM_READ_DIM(strlen)
        M_ALLOCATE(character(len=strlen) :: values(n_pairs))

        ! Read each line and parse key-value pairs
        do i_pair = 1, n_pairs
            CM_READ(keys(i_pair))
            CM_READ(values(i_pair))
        end do

        close (unit)
        if (DEBUG) print *, "Read manifest with ", n_pairs, " entries"
    end subroutine read_manifest_generic

    !> AUTHOR_AARON_SCHROEDER
    !| Extract and parse the manifest file for generic key-value pairs.
    subroutine extract_and_parse_manifest(zip_handle, keys, filenames, ierr)
        use tox_conversions, only: int32_as_c_size
        type(c_ptr), intent(in) :: zip_handle
            !! Zip file connection
        character(len=:), allocatable, intent(out) :: keys(:)
            !! Array of keys from manifest
        character(len=:), allocatable, intent(out) :: filenames(:)
            !! Array of filenames from manifest
        integer(int32), intent(out) :: ierr
            !! Error code

        type(c_ptr) :: file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: bytes_read
        integer(int32) :: unit, iostat
        integer(int32), parameter :: CHUNK_SIZE = 4096
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        integer(c_size_t) :: chunk_size_c

        ! Extract and parse the manifest file
        call extract_file_from_zip(zip_handle, "manifest.txt", ierr)
        if (is_err(ierr)) return
        call read_manifest_generic("manifest.txt", keys, filenames, ierr)

        if (is_err(ierr)) then
            if (DEBUG) print *, "Error parsing manifest file"
            return
        end if
    end subroutine extract_and_parse_manifest

    !> AUTHOR_AARON_SCHROEDER
    !| Save standard tox data
    subroutine save_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, &
                             expression_file, gene_to_family, gene_to_family_file, &
                             family_ids, family_ids_file, family_centroids, &
                             family_centroids_file, shift_vectors, shift_vectors_file)
        implicit none

        character(len=*), intent(in) :: zip_filename
        !! Zip filename
        character(len=*), intent(in), optional :: gene_ids(:)
        !! Gene ids array, will be saved if provided
        character(len=*), intent(in), optional :: family_ids(:)
        !! Family ids array, will be saved if provided
        real(real64), intent(in), optional :: expression(:, :)
        !! Expression vectors array, will be saved if provided
        real(real64), intent(in), optional :: family_centroids(:, :)
        !! Family centroids array, will be saved if provided
        real(real64), intent(in), optional :: shift_vectors(:, :)
        !! Shift vectors array, will be saved if provided
        integer(int32), intent(in), optional :: gene_to_family(:)
        !! Gene to family mapping array, will be saved if provided
        character(len=*), intent(in), optional :: gene_ids_file
        !! Name of the gene ids file
        character(len=*), intent(in), optional :: expression_file
        !! Name of the expression file
        character(len=*), intent(in), optional :: gene_to_family_file
        !! Name of the gene to family mapping file
        character(len=*), intent(in), optional :: family_ids_file
        !! Name of the family ids file
        character(len=*), intent(in), optional :: family_centroids_file
        !! Name of the family centroids file
        character(len=*), intent(in), optional :: shift_vectors_file
        !! Name of the shift vectors file
        integer(int32), intent(out) :: ierr
        !! Error code

        character(len=:), allocatable :: actual_gene_ids_file, actual_expression_file, actual_gene_to_family_file, &
                                         actual_family_ids_file, actual_family_centroids_file, actual_shift_vectors_file
        logical :: gene_ids_present, expression_present, gene_to_family_present, &
                   family_ids_present, family_centroids_present, shift_vectors_present
        integer(int32) :: temp_ierr
        character(len=:), allocatable :: keys(:), filenames(:)
        integer :: count, i

        call set_ok(ierr)
        call set_ok(temp_ierr)

        ! Determine which arrays are present
        gene_ids_present = present(gene_ids) .and. present(gene_ids_file)
        expression_present = present(expression) .and. present(expression_file)
        gene_to_family_present = present(gene_to_family) .and. present(gene_to_family_file)
        family_ids_present = present(family_ids) .and. present(family_ids_file)
        family_centroids_present = present(family_centroids) .and. present(family_centroids_file)
        shift_vectors_present = present(shift_vectors) .and. present(shift_vectors_file)

        ! Count the number of present arrays to determine the size of keys and filenames
        count = 0
        if (gene_ids_present) count = count + 1
        if (expression_present) count = count + 1
        if (gene_to_family_present) count = count + 1
        if (family_ids_present) count = count + 1
        if (family_centroids_present) count = count + 1
        if (shift_vectors_present) count = count + 1

        ! Allocate keys and filenames arrays
        M_ALLOCATE(character(len=32) :: keys(count))
        M_ALLOCATE(character(len=256) :: filenames(count))

        ! Save data files and populate keys and filenames
        i = 1
        if (gene_ids_present) then
            actual_gene_ids_file = gene_ids_file
            call save_gene_ids(gene_ids, actual_gene_ids_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'gene_ids'
                filenames(i) = actual_gene_ids_file
                i = i + 1
            else
                return
            end if
        else
            actual_gene_ids_file = ""
        end if

        if (expression_present) then
            actual_expression_file = expression_file
            call save_expression_vectors(expression, actual_expression_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'expression'
                filenames(i) = actual_expression_file
                i = i + 1
            else
                return
            end if
        else
            actual_expression_file = ""
        end if

        if (gene_to_family_present) then
            actual_gene_to_family_file = gene_to_family_file
            call save_gene_to_family(gene_to_family, actual_gene_to_family_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'gene_to_family'
                filenames(i) = actual_gene_to_family_file
                i = i + 1
            else
                return
            end if
        else
            actual_gene_to_family_file = ""
        end if

        if (family_ids_present) then
            actual_family_ids_file = family_ids_file
            call save_family_ids(family_ids, actual_family_ids_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'family_ids'
                filenames(i) = actual_family_ids_file
                i = i + 1
            else
                return
            end if
        else
            actual_family_ids_file = ""
        end if

        if (family_centroids_present) then
            actual_family_centroids_file = family_centroids_file
            call save_family_centroids(family_centroids, actual_family_centroids_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'family_centroids'
                filenames(i) = actual_family_centroids_file
                i = i + 1
            else
                return
            end if
        else
            actual_family_centroids_file = ""
        end if

        if (shift_vectors_present) then
            actual_shift_vectors_file = shift_vectors_file
            call save_shift_vectors(shift_vectors, actual_shift_vectors_file, ierr)
            if (is_ok(ierr)) then
                keys(i) = 'shift_vectors'
                filenames(i) = actual_shift_vectors_file
                i = i + 1
            else
                return
            end if
        else
            actual_shift_vectors_file = ""
        end if

        if (is_ok(ierr)) then
            call create_zip_archive(zip_filename, keys, filenames, ierr)
        end if

        ! Clean up temporary files
        call cleanup_temporary_files(gene_ids_present, actual_gene_ids_file, "Gene IDs")
        call cleanup_temporary_files(expression_present, actual_expression_file, "Expression")
        call cleanup_temporary_files(gene_to_family_present, actual_gene_to_family_file, "Gene to family mapping")
        call cleanup_temporary_files(family_ids_present, actual_family_ids_file, "Family IDs")
        call cleanup_temporary_files(family_centroids_present, actual_family_centroids_file, "Centroids")
        call cleanup_temporary_files(shift_vectors_present, actual_shift_vectors_file, "Shift vectors")

        deallocate (keys, filenames)

    contains
        !| Deletes one intermediate array file after it has been added to the archive. Failure to
        !| remove it is only ever a warning (via the enclosing subroutine's `temp_ierr`, never
        !| `ierr`), since a leftover temp file does not affect whether the archive itself was
        !| written successfully.
        subroutine cleanup_temporary_files(file_present, filename, description)
            logical, intent(in) :: file_present
                !! Whether this array was actually saved (and thus has a temp file to remove)
            character(len=*), intent(in) :: filename
                !! Path of the temporary file to delete
            character(len=*), intent(in) :: description
                !! Human-readable label used in the warning message if deletion fails

            if (file_present .and. len_trim(filename) > 0) then
                call delete_file(filename, temp_ierr)
                if (is_err(temp_ierr)) then
                    if (DEBUG) write (*, *) 'Warning: ', trim(description), ' file could not be removed: ', trim(filename)
                end if
            end if
        end subroutine cleanup_temporary_files
    end subroutine save_tox_data

    !> AUTHOR_AARON_SCHROEDER
    !| Read standard tox from a zip archive
    subroutine read_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, expression_file, &
                             gene_to_family, gene_to_family_file, family_ids, family_ids_file, &
                             family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
        use f42_serde_arrays_utils, only: get_array_metadata, REAL_TYPE_CODE, INTEGER_TYPE_CODE
        use tox_data_read_write
        use iso_fortran_env, only: real64, int32
        implicit none

        character(len=*), intent(in) :: zip_filename
        !! Name of the zipfile
        integer(int32), intent(out) :: ierr
        !! Error code
        character(len=:), allocatable, optional, intent(out) :: gene_ids(:)
        !! Gene IDs array, will be populated if provided
        character(len=:), allocatable, optional, intent(out) :: family_ids(:)
        !! Family IDs array, will be populated if provided
        real(real64), allocatable, optional, intent(out) :: expression(:, :)
        !! Expression vectors array, will be populated if provided
        real(real64), allocatable, optional, intent(out) :: family_centroids(:, :)
        !! Family centroids array, will be populated if provided
        real(real64), allocatable, optional, intent(out) :: shift_vectors(:, :)
        !! Shift vectors array, will be populated if provided
        integer(int32), allocatable, optional, intent(out) :: gene_to_family(:)
        !! Gene to family mapping array, will be populated if provided
        character(len=:), allocatable, optional, intent(out) :: gene_ids_file
        !! Name of the gene ids file in the zip archive
        character(len=:), allocatable, optional, intent(out) :: expression_file
        !! Name of the expression vectors file in the zip archive
        character(len=:), allocatable, optional, intent(out) :: gene_to_family_file
        !! Name of the gene to family mapping file in the zip archive
        character(len=:), allocatable, optional, intent(out) :: family_ids_file
        !! Name of the family ids file in the zip archive
        character(len=:), allocatable, optional, intent(out) :: family_centroids_file
        !! Name of the family centroids file in the zip archive
        character(len=:), allocatable, optional, intent(out) :: shift_vectors_file
        !! Name of the shift vectors file in the zip archive

        character(len=:), allocatable :: keys(:), filenames(:)
        integer(int32) :: i_key, type_code
        logical :: gene_ids_requested, expression_requested, gene_to_family_requested, &
                   family_ids_requested, family_centroids_requested, shift_vectors_requested
        integer(int32) :: max_dims, ndims, dims(5)
        character(len=:), allocatable :: extracted_gene_ids_file, extracted_expression_file, &
                                         extracted_gene_to_family_file, extracted_family_ids_file, &
                                         extracted_family_centroids_file, extracted_shift_vectors_file
        integer(int32) :: temp_ierr
        ! Sanity bounds for dims/type_code read from on-disk archive member metadata, to guard
        ! against negative or absurdly large values from a corrupted/malicious archive
        integer(int32), parameter :: MAX_REASONABLE_DIM = 100000000
        integer(int32), parameter :: MAX_REASONABLE_CHARLEN = 100000

        call set_ok(ierr)
        max_dims = 5

        if (DEBUG) write (*, *) 'Extracting zip archive...'

        call extract_zip_archive(zip_filename, keys, filenames, ierr)
        if (is_err(ierr)) return

        ! Find standard files by their keys
        extracted_gene_ids_file = ""
        extracted_expression_file = ""
        extracted_gene_to_family_file = ""
        extracted_family_ids_file = ""
        extracted_family_centroids_file = ""
        extracted_shift_vectors_file = ""

        ! This key vocabulary is the other half of the contract written by
        ! [[tox_data_archive(module):save_tox_data(subroutine)]] -- the two must be kept in sync,
        ! since a key renamed on one side but not the other would silently fall through to the
        ! "non-standard key" branch below instead of failing loudly.
        do i_key = 1, size(keys)
            select case (trim(keys(i_key)))
            case ('gene_ids')
                extracted_gene_ids_file = trim(filenames(i_key))
            case ('expression')
                extracted_expression_file = trim(filenames(i_key))
            case ('gene_to_family')
                extracted_gene_to_family_file = trim(filenames(i_key))
            case ('family_ids')
                extracted_family_ids_file = trim(filenames(i_key))
            case ('family_centroids')
                extracted_family_centroids_file = trim(filenames(i_key))
            case ('shift_vectors')
                extracted_shift_vectors_file = trim(filenames(i_key))
            case default
                print *, "Found non-standard key in archive: ", trim(keys(i_key)), " in file: ", trim(filenames(i_key))
            end select
        end do

        ! Return filenames if requested
        if (present(gene_ids_file)) gene_ids_file = extracted_gene_ids_file
        if (present(expression_file)) expression_file = extracted_expression_file
        if (present(gene_to_family_file)) gene_to_family_file = extracted_gene_to_family_file
        if (present(family_ids_file)) family_ids_file = extracted_family_ids_file
        if (present(family_centroids_file)) family_centroids_file = extracted_family_centroids_file
        if (present(shift_vectors_file)) shift_vectors_file = extracted_shift_vectors_file

        ! Load the arrays that are requested and available
        !
        ! NOTE on `type_code` below: for numeric arrays it is one of the negative
        ! M_INTEGER_TYPE_CODE / M_REAL_TYPE_CODE sentinels (see f42_serde_arrays_utils), but for
        ! character arrays the on-disk header instead stores the fixed string length in that same
        ! field (see serialize_char_helper), which is always >= 0. That is why the character
        ! branches below check `type_code >= 0` (treating it as a length) while the numeric
        ! branches check `type_code == REAL_TYPE_CODE`/`INTEGER_TYPE_CODE` (treating it as a
        ! sentinel) -- the two checks are not testing the same kind of value.
        gene_ids_requested = present(gene_ids) .and. len_trim(extracted_gene_ids_file) > 0
        if (gene_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_gene_ids_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 1 .and. type_code >= 0 .and. type_code <= MAX_REASONABLE_CHARLEN &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata with proper character length
                M_ALLOCATE(character(len=type_code) :: gene_ids(dims(1)))
                call load_gene_ids(gene_ids, extracted_gene_ids_file, ierr)
                if (is_err(ierr)) return
            else
                print *, "Error getting metadata for gene_ids file"
                return
            end if
        end if

        expression_requested = present(expression) .and. len_trim(extracted_expression_file) > 0
        if (expression_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_expression_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 2 .and. type_code == REAL_TYPE_CODE &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM &
                .and. dims(2) >= 0 .and. dims(2) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata
                M_ALLOCATE(expression(dims(1), dims(2)))
                call load_expression_vectors(expression, extracted_expression_file, ierr)
                if (is_err(ierr)) return
            else
                print *, "Error getting metadata for expression file"
                return
            end if
        end if

        gene_to_family_requested = present(gene_to_family) .and. len_trim(extracted_gene_to_family_file) > 0
        if (gene_to_family_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_gene_to_family_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 1 .and. type_code == INTEGER_TYPE_CODE &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata
                M_ALLOCATE(gene_to_family(dims(1)))
                call load_gene_to_family(gene_to_family, extracted_gene_to_family_file, ierr)
                if (is_err(ierr)) return
            else
                print *, "Error getting metadata for gene_to_family file"
                return
            end if
        end if

        family_ids_requested = present(family_ids) .and. len_trim(extracted_family_ids_file) > 0
        if (family_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_family_ids_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 1 .and. type_code >= 0 .and. type_code <= MAX_REASONABLE_CHARLEN &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata with proper character length
                M_ALLOCATE(character(len=type_code) :: family_ids(dims(1)))
                call load_family_ids(family_ids, extracted_family_ids_file, ierr)
            else
                print *, "Error getting metadata for family_ids file"
                return
            end if
        end if

        family_centroids_requested = present(family_centroids) .and. len_trim(extracted_family_centroids_file) > 0
        if (family_centroids_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_family_centroids_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 2 .and. type_code == REAL_TYPE_CODE &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM &
                .and. dims(2) >= 0 .and. dims(2) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata
                M_ALLOCATE(family_centroids(dims(1), dims(2)))
                call load_family_centroids(family_centroids, extracted_family_centroids_file, ierr)
                if (is_err(ierr)) return
            else
                print *, "Error getting metadata for family_centroids file"
                return
            end if
        end if

        shift_vectors_requested = present(shift_vectors) .and. len_trim(extracted_shift_vectors_file) > 0
        if (shift_vectors_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_shift_vectors_file, dims, max_dims, ndims, type_code, ierr)
            if (is_ok(ierr) .and. ndims == 2 .and. type_code == REAL_TYPE_CODE &
                .and. dims(1) >= 0 .and. dims(1) <= MAX_REASONABLE_DIM &
                .and. dims(2) >= 0 .and. dims(2) <= MAX_REASONABLE_DIM) then
                ! Allocate array based on metadata
                M_ALLOCATE(shift_vectors(dims(1), dims(2)))
                call load_shift_vectors(shift_vectors, extracted_shift_vectors_file, ierr)
                if (is_err(ierr)) return
            else
                print *, "Error getting metadata for shift_vectors file"
                return
            end if
        end if

        ! Use a separate status variable for cleanup so a temp-file deletion failure does not
        ! mask the overall success of the actual data load (mirrors save_tox_data's temp_ierr)
        call set_ok(temp_ierr)
        if (len_trim(extracted_gene_ids_file) > 0) call delete_file(extracted_gene_ids_file, temp_ierr)
        if (len_trim(extracted_expression_file) > 0) call delete_file(extracted_expression_file, temp_ierr)
        if (len_trim(extracted_gene_to_family_file) > 0) call delete_file(extracted_gene_to_family_file, temp_ierr)
        if (len_trim(extracted_family_ids_file) > 0) call delete_file(extracted_family_ids_file, temp_ierr)
        if (len_trim(extracted_family_centroids_file) > 0) call delete_file(extracted_family_centroids_file, temp_ierr)
        if (len_trim(extracted_shift_vectors_file) > 0) call delete_file(extracted_shift_vectors_file, temp_ierr)
        call delete_file("manifest.txt", temp_ierr)
        if (DEBUG .and. is_err(temp_ierr)) print *, "Warning: failed to clean up one or more temporary extracted files"

        deallocate (keys, filenames)
    end subroutine read_tox_data

end module tox_data_archive

!> C binding for generic archive creation with arrays of keys and filenames
subroutine create_zip_archive_c(zip_filename, zip_len, &
                                keys, keys_len, keys_count, &
                                filenames, filenames_len, filenames_count, &
                                ierr) bind(C, name="create_zip_archive_c")
    use tox_data_archive, only: create_zip_archive
    use tox_conversions, only: c_char_2d_as_string, c_char_1d_as_string
    use iso_c_binding, only: c_int, c_char
    use tox_errors, only: is_err, set_ok, set_err_once, ERR_INVALID_INPUT
    use iso_fortran_env, only: int32
    M_USE_NULL_VALIDATION

    ! Input arguments
    integer(c_int), intent(in), target :: zip_len
    !! Length of the zip filename
    character(kind=c_char, len=1), intent(in), target :: zip_filename(zip_len)
    !! Zip Filename as c_chars
    integer(c_int), intent(in), target :: keys_count
    !! number of keys
    integer(c_int), intent(in), target :: keys_len
    !! lengths of the keys
    character(kind=c_char, len=1), intent(in), target :: keys(keys_len, keys_count)
    !! Keys as c_chars
    integer(c_int), intent(in), target :: filenames_count
    !! Number of files
    integer(c_int), intent(in), target :: filenames_len
    !! Length of the filenames
    character(kind=c_char, len=1), intent(in), target :: filenames(filenames_len, filenames_count)
    !! Filenames as c_chars
    integer(c_int), intent(out), target :: ierr
    !! Error code

    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: f_keys(:)
    character(len=:), allocatable :: f_filenames(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(zip_len)
    M_CHECK_NON_NULL(keys_count)
    M_CHECK_NON_NULL(keys_len)
    M_CHECK_NON_NULL(filenames_count)
    M_CHECK_NON_NULL(filenames_len)
    M_CHECK_NON_NULL(zip_filename)
    M_CHECK_NON_NULL(keys)
    M_CHECK_NON_NULL(filenames)

    call set_ok(ierr)

    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if (is_err(ierr)) return

    call c_char_2d_as_string(keys, f_keys, ierr)
    if (is_err(ierr)) return

    call c_char_2d_as_string(filenames, f_filenames, ierr)
    if (is_err(ierr)) return

    ! Validate array sizes
    if (size(f_keys) /= size(f_filenames)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    call create_zip_archive(f_zip_filename, f_keys, f_filenames, ierr)

end subroutine create_zip_archive_c

!> C binding for extract_zip_archive - can be called directly from Python via ctypes
subroutine extract_zip_archive_c(zip_filename, filename_len, ierr) &
    bind(C, name="extract_zip_archive_c")
    use tox_data_archive, only: extract_zip_archive
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, is_err
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32
    M_USE_NULL_VALIDATION

    ! Input arguments
    integer(c_int), intent(in), target :: filename_len
    !! Length of the filename
    character(kind=c_char, len=1), intent(in), target :: zip_filename(filename_len)
    !! Zip filename length
    integer(c_int), intent(out), target :: ierr
    !! Error code

    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: keys(:), filenames(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(filename_len)
    M_CHECK_NON_NULL(zip_filename)

    call set_ok(ierr)

    ! Convert C string to Fortran string
    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if (is_err(ierr)) return

    call extract_zip_archive(f_zip_filename, keys, filenames, ierr)

end subroutine extract_zip_archive_c
