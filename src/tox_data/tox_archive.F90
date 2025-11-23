module tox_archive
    use safeguard
    use iso_c_binding, only: c_ptr, c_char, c_int, c_int64_t, c_size_t, c_signed_char, c_f_pointer, c_loc, c_associated, c_null_char, c_null_ptr
    use tox_data_read_write
    use tox_errors, only: set_ok, set_err_once, is_err, ERR_FILE_OPEN, ERR_ALLOC_FAIL, ERR_FILE_ADD
    use tox_errors, only: ERR_FILE_CLOSE, ERR_FILE_EXTRACT, ERR_INVALID_INPUT
    use tox_errors, only: ERR_POINTER_NULL, ERR_WRITE_DATA, ERR_READ_DATA, ERR_MISSING_MANIFEST
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
    integer(int32), parameter :: DATA_TYPE_FILE   = 1
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

    !> Creates a zip archive with generic file lists
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
            if(DEBUG) print *, "Error: keys and filenames arrays must have same size"
            return
        end if

        if (len_trim(zip_filename) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if(DEBUG) print *, "Error: zip_filename cannot be empty"
            return
        end if
        
        ! Open ZIP archive
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
        if(is_err(ierr)) then
            error = zip_close(zip_handle)
            return
        end if
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (is_err(error)) then
            call set_err_once(ierr, error)
            if(DEBUG) print *, "Error closing ZIP file: ", error
        else if (is_ok(ierr)) then
            if(DEBUG) print *, "ZIP archive created successfully: ", trim(zip_filename)
        end if
    end subroutine create_zip_archive

    !> Extract a zip archive and return all key-value pairs from manifest
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
        allocate(character(len=256) :: keys(0), filenames(0))

        call set_ok(ierr)
        call set_ok(error)

        ! Check if file exists
        inquire(file=zip_filename, exist=file_exists)
        if (.not. file_exists) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if(DEBUG) print *, "ZIP file does not exist: ", trim(zip_filename)
            return
        end if

        ! Open ZIP archive
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_READ_ONLY, error)
        if (error /= 0 .or. .not. c_associated(zip_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if(DEBUG) print *, "Error opening ZIP file for reading: ", error
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
            if(DEBUG) print *, "Error closing ZIP file: ", error
            call set_err_once(ierr, ERR_FILE_CLOSE)
            return
        end if

        if(DEBUG) print *, "ZIP archive extracted successfully: ", trim(zip_filename)
    end subroutine extract_zip_archive

    !> Delete a file from the disk
    subroutine delete_file(filename, ierr)
        character(len=*), intent(in) :: filename
            !! File to delete
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32) :: unit, iostat

        call set_ok(iostat)
        call set_ok(ierr)
        
        open(newunit=unit, file=filename, iostat=iostat, status='old')
        if (is_ok(iostat)) then
            close(unit, status='delete')
        else 
            call set_err_once(ierr, ERR_FILE_OPEN)
        end if
    end subroutine delete_file

    !> Helper function to get the name of a ZIP entry
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
        ! Get name from ZIP
        call c_f_pointer(zip_get_name(zip_handle, entry_index, 0), name_ptr, [MAX_NAME_LENGTH])
        call c_char_1d_as_string(name_ptr, entry_name, ierr)
    end subroutine get_zip_entry_name

    ! Unified subroutine to extract a file from ZIP archive
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
            if(DEBUG) print *, "Error opening file in ZIP: ", trim(filename)
            return
        end if
        
        ! Open output file
        open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
            iostat=iostat, status='replace', action='write')
        if (is_err(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if(DEBUG) print *, "Error creating file: ", trim(filename)
            error = zip_fclose(file_handle)
            return
        end if
        
        ! Read and write in chunks
        allocate(buffer(CHUNK_SIZE), stat=iostat)
        if (is_err(iostat)) then
            call set_err_once(ierr, ERR_ALLOC_FAIL)
            if(DEBUG) print *, "Error allocating buffer for: ", trim(filename)
            close(unit, status='delete')
            error = zip_fclose(file_handle)
            return
        end if

        call int32_as_c_size(CHUNK_SIZE, chunk_size_c)
        do
            bytes_read = zip_fread(file_handle, c_loc(buffer), chunk_size_c)
            if (bytes_read <= 0) exit
            if (bytes_read > 0) then
                write(unit, iostat=iostat) buffer(1:bytes_read)
                if (is_err(iostat)) then
                    call set_err_once(ierr, ERR_WRITE_DATA)
                    if(DEBUG) print *, "Error writing file: ", trim(filename)
                    exit
                end if
            end if
        end do
        
        ! Clean up
        if (allocated(buffer)) deallocate(buffer)
        close(unit)
        error = zip_fclose(file_handle)
        
        if (is_err(error)) then
            call set_err_once(ierr, ERR_FILE_CLOSE)
            if(DEBUG) print *, "Error closing file in ZIP: ", trim(filename)
        end if
        
        if (is_ok(ierr)) then
            if(DEBUG) print *, "Extracted: ", trim(filename)
        end if
    end subroutine extract_file_from_zip

    ! Unified subroutine to add data to ZIP (handles both files and strings)
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

        if (len_trim(data_source) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        select case (data_type)

        case (DATA_TYPE_FILE)
            ! Open file
            open(newunit=unit, file=data_source, access='stream', form='unformatted', &
                iostat=iostat, status='old')
            if (is_err(iostat)) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                if (DEBUG) print *, "Error opening file: ", trim(data_source)
                return
            end if

            inquire(unit, size=file_size)
            call int32_as_c_size(file_size, data_len)

            if (file_size == 0) then
                close(unit)
                source = zip_source_buffer(zip_handle, c_null_ptr, 0_c_size_t, 0)
            else
                c_data = malloc(data_len)
                if (.not. c_associated(c_data)) then
                    call set_err_once(ierr, ERR_POINTER_NULL)
                    close(unit)
                    return
                end if

                call c_f_pointer(c_data, file_data, [file_size])
                read(unit, iostat=iostat) file_data
                close(unit)

                if (is_err(iostat)) then
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
            if (data_type == DATA_TYPE_FILE .or. data_type == DATA_TYPE_STRING) call free(c_data)
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

    !> Write manifest from given key-value pairs
    subroutine write_manifest(keys, filenames, manifest_filename, ierr)
        character(len=*), intent(in) :: keys(:)
            !! Array of keys for manifest entries
        character(len=*), intent(in) :: filenames(:)
            !! Array of filenames for manifest entries  
        character(len=*), intent(in) :: manifest_filename
            !! Name of the manifest file (should be manifest.txt)
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(int32) :: unit, iostat, i
        character(len=:), allocatable :: line

        call set_ok(ierr)
        
        ! Validate input arrays have same size
        if (size(keys) /= size(filenames)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if(DEBUG) print *, "Error: keys and filenames arrays must have same size"
            return
        end if
        
        ! Open the manifest file for writing
        open(newunit=unit, file=manifest_filename, status='replace', iostat=iostat)
        if (is_err(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if(DEBUG) print *, "Error creating manifest file: ", trim(manifest_filename)
            return
        end if
        
        ! Write each key-value pair to the manifest
        do i = 1, size(keys)
            if (len_trim(keys(i)) > 0 .and. len_trim(filenames(i)) > 0) then
                line = trim(keys(i)) // '=' // trim(filenames(i))
                write(unit, '(a)') trim(line)
            end if
        end do
        
        ! Close the manifest file
        close(unit)
        
        if(DEBUG) print *, "Manifest created successfully with ", size(keys), " entries"
    end subroutine write_manifest

    !> Read manifest file and return key-value pairs
    subroutine read_manifest_generic(manifest_filename, keys, values, ierr)
        character(len=*), intent(in) :: manifest_filename
            !! Filename of the manifest (should be manifest.txt)
        character(len=:), allocatable, intent(out) :: keys(:)
            !! Array of keys from manifest
        character(len=:), allocatable, intent(out) :: values(:)
            !! Array of values from manifest
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(int32) :: unit, iostat, line_count, i, eq_pos
        character(len=1024) :: line  ! Increased buffer size
        character(len=:), allocatable :: temp_keys(:), temp_values(:)
        integer(int32), parameter :: MAX_LINES = 100  ! Reasonable maximum
        
        call set_ok(ierr)
        
        ! Initialize temporary arrays with proper length
        allocate(character(len=256) :: temp_keys(MAX_LINES))
        allocate(character(len=256) :: temp_values(MAX_LINES))
        line_count = 0
        
        ! Open the manifest file for reading
        open(newunit=unit, file=manifest_filename, status='old', iostat=iostat, action='read')
        if (is_err(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if(DEBUG) print *, "Error opening manifest file: ", trim(manifest_filename)
            deallocate(temp_keys, temp_values)
            return
        end if
        
        ! Read each line and parse key-value pairs
        do i = 1, MAX_LINES
            read(unit, '(a)', iostat=iostat) line
            if (iostat == iostat_end) exit  ! End of file
            if (is_err(iostat)) then
                call set_err_once(ierr, ERR_READ_DATA)
                if(DEBUG) print *, "Error reading manifest file: ", trim(manifest_filename)
                exit
            end if
            
            ! Skip empty lines and lines without '='
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            
            eq_pos = index(line, '=')
            if (eq_pos == 0) cycle  ! Skip lines without '='
            
            ! Extract key and value
            temp_keys(line_count + 1) = trim(adjustl(line(1:eq_pos-1)))
            temp_values(line_count + 1) = trim(adjustl(line(eq_pos+1:)))
            
            ! Only count if both key and value are non-empty
            if (len_trim(temp_keys(line_count + 1)) > 0 .and. &
                len_trim(temp_values(line_count + 1)) > 0) then
                line_count = line_count + 1
            end if
        end do
        
        close(unit)
        
        ! Allocate output arrays with correct size
        if (line_count > 0) then
            allocate(character(len=256) :: keys(line_count), values(line_count))
            do i = 1, line_count
                keys(i) = trim(temp_keys(i))
                values(i) = trim(temp_values(i))
            end do
        else
            allocate(character(len=256) :: keys(0), values(0))
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if(DEBUG) print *, "No valid key-value pairs found in manifest"
        end if
        
        deallocate(temp_keys, temp_values)
        
        if(DEBUG) print *, "Read manifest with ", line_count, " entries"
    end subroutine read_manifest_generic

    !> Extract and parse the manifest file for generic key-value pairs.
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
        file_handle = zip_fopen(zip_handle, "manifest.txt"//c_null_char, 0)
        if (c_associated(file_handle)) then
            ! Open output file for manifest
            open(newunit=unit, file="manifest.txt", access='stream', form='unformatted', &
                iostat=iostat, status='replace', action='write')
            if (is_err(iostat)) then
                if(DEBUG) print *, "Error creating manifest file"
                error = zip_fclose(file_handle)
                call set_err_once(ierr, ERR_FILE_OPEN)
                return
            end if

            ! Read and write in chunks
            allocate(buffer(CHUNK_SIZE), stat=iostat)
            if (is_ok(iostat)) then
                call int32_as_c_size(CHUNK_SIZE, chunk_size_c)
                do
                    bytes_read = zip_fread(file_handle, c_loc(buffer), chunk_size_c)
                    if (bytes_read <= 0) exit
                    write(unit, iostat=iostat) buffer(1:bytes_read)
                    if (is_err(iostat)) then
                        if(DEBUG) print *, "Error writing manifest file"
                        call set_err_once(ierr, ERR_WRITE_DATA)
                        exit
                    end if
                end do
            else
                call set_err_once(ierr, ERR_ALLOC_FAIL)
                error = zip_close(zip_handle)
                return
            end if

            ! Clean up manifest extraction
            if (allocated(buffer)) deallocate(buffer)
            close(unit)
            error = zip_fclose(file_handle)
            if(is_err(error)) then
                if(DEBUG) print *, "Error closing manifest file in ZIP"
                call set_err_once(ierr, ERR_FILE_CLOSE)
                return
            end if

            ! Parse the manifest file
            call read_manifest_generic("manifest.txt", keys, filenames, ierr)

            if (is_err(ierr)) then
                if(DEBUG) print *, "Error parsing manifest file"
                return
            end if

        else
            if(DEBUG) print *, "No manifest file found in ZIP archive"
            call set_err_once(ierr, ERR_MISSING_MANIFEST)
        end if
    end subroutine extract_and_parse_manifest

    !> Save standard tox data
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
        real(real64), intent(in), optional :: expression(:,:)
        !! Expression vectors array, will be saved if provided
        real(real64), intent(in), optional :: family_centroids(:,:)
        !! Family centroids array, will be saved if provided
        real(real64), intent(in), optional :: shift_vectors(:,:)
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
        allocate(character(len=32) :: keys(count))
        allocate(character(len=256) :: filenames(count))

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

        deallocate(keys, filenames)

    contains
        subroutine cleanup_temporary_files(file_present, filename, description)
            logical, intent(in) :: file_present
            character(len=*), intent(in) :: filename
            character(len=*), intent(in) :: description

            if (file_present .and. len_trim(filename) > 0) then
                call delete_file(filename, temp_ierr)
                if(is_err(temp_ierr)) then
                    if(DEBUG) write(*,*) 'Warning: ', trim(description), ' file could not be removed: ', trim(filename)
                end if
            end if
        end subroutine cleanup_temporary_files
    end subroutine save_tox_data

    !> Read standard tox from a zip archive
    subroutine read_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, expression_file, &
                        gene_to_family, gene_to_family_file, family_ids, family_ids_file, &
                        family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
        use array_utils, only: get_array_metadata
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
        real(real64), allocatable, optional, intent(out) :: expression(:,:)
        !! Expression vectors array, will be populated if provided
        real(real64), allocatable, optional, intent(out) :: family_centroids(:,:)
        !! Family centroids array, will be populated if provided
        real(real64), allocatable, optional, intent(out) :: shift_vectors(:,:)
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
        integer(int32) :: i
        logical :: gene_ids_requested, expression_requested, gene_to_family_requested, &
                family_ids_requested, family_centroids_requested, shift_vectors_requested
        integer(int32) :: max_dims, ndims, dims(5), char_len
        character(len=:), allocatable :: extracted_gene_ids_file, extracted_expression_file, &
                                    extracted_gene_to_family_file, extracted_family_ids_file, &
                                    extracted_family_centroids_file, extracted_shift_vectors_file
        
        call set_ok(ierr)
        max_dims = 5
        
        if(DEBUG) write(*,*) 'Extracting zip archive...'

        call extract_zip_archive(zip_filename, keys, filenames, ierr)
        if (is_err(ierr)) return
        
        ! Find standard files by their keys
        extracted_gene_ids_file = ""
        extracted_expression_file = ""
        extracted_gene_to_family_file = ""
        extracted_family_ids_file = ""
        extracted_family_centroids_file = ""
        extracted_shift_vectors_file = ""
        
        do i = 1, size(keys)
            select case (trim(keys(i)))
                case ('gene_ids')
                    extracted_gene_ids_file = trim(filenames(i))
                case ('expression')
                    extracted_expression_file = trim(filenames(i))
                case ('gene_to_family')
                    extracted_gene_to_family_file = trim(filenames(i))
                case ('family_ids')
                    extracted_family_ids_file = trim(filenames(i))
                case ('family_centroids')
                    extracted_family_centroids_file = trim(filenames(i))
                case ('shift_vectors')
                    extracted_shift_vectors_file = trim(filenames(i))
                case default
                    print *, "Found non-standard key in archive: ", trim(keys(i)), " in file: ", trim(filenames(i))
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
        gene_ids_requested = present(gene_ids) .and. len_trim(extracted_gene_ids_file) > 0
        if (gene_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_gene_ids_file, dims, max_dims, ndims, ierr, char_len)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: gene_ids(dims(1)))
                call load_gene_ids(gene_ids, extracted_gene_ids_file, ierr)
                if(is_err(ierr)) return
            else
                print *, "Error getting metadata for gene_ids file"
                return
            end if
        end if
        
        expression_requested = present(expression) .and. len_trim(extracted_expression_file) > 0
        if (expression_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_expression_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(expression(dims(1), dims(2)))
                call load_expression_vectors(expression, extracted_expression_file, ierr)
                if(is_err(ierr)) return
            else
                print *, "Error getting metadata for expression file"
                return
            end if
        end if
        
        gene_to_family_requested = present(gene_to_family) .and. len_trim(extracted_gene_to_family_file) > 0
        if (gene_to_family_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_gene_to_family_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata
                allocate(gene_to_family(dims(1)))
                call load_gene_to_family(gene_to_family, extracted_gene_to_family_file, ierr)
                if(is_err(ierr)) return
            else
                print *, "Error getting metadata for gene_to_family file"
                return
            end if
        end if
        
        family_ids_requested = present(family_ids) .and. len_trim(extracted_family_ids_file) > 0
        if (family_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_family_ids_file, dims, max_dims, ndims, ierr, char_len)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: family_ids(dims(1)))
                call load_family_ids(family_ids, extracted_family_ids_file, ierr)
            else
                print *, "Error getting metadata for family_ids file"
                return
            end if
        end if
        
        family_centroids_requested = present(family_centroids) .and. len_trim(extracted_family_centroids_file) > 0
        if (family_centroids_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_family_centroids_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(family_centroids(dims(1), dims(2)))
                call load_family_centroids(family_centroids, extracted_family_centroids_file, ierr)
                if(is_err(ierr)) return
            else
                print *, "Error getting metadata for family_centroids file"
                return
            end if
        end if
        
        shift_vectors_requested = present(shift_vectors) .and. len_trim(extracted_shift_vectors_file) > 0
        if (shift_vectors_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_shift_vectors_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr).and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(shift_vectors(dims(1), dims(2)))
                call load_shift_vectors(shift_vectors, extracted_shift_vectors_file, ierr)
                if(is_err(ierr)) return
            else
                print *, "Error getting metadata for shift_vectors file"
                return
            end if
        end if

        if(len_trim(extracted_gene_ids_file) > 0) call delete_file(extracted_gene_ids_file, ierr)
        if(len_trim(extracted_expression_file) > 0) call delete_file(extracted_expression_file, ierr)
        if(len_trim(extracted_gene_to_family_file) > 0) call delete_file(extracted_gene_to_family_file, ierr)
        if(len_trim(extracted_family_ids_file) > 0) call delete_file(extracted_family_ids_file, ierr)
        if(len_trim(extracted_family_centroids_file) > 0) call delete_file(extracted_family_centroids_file, ierr)
        if(len_trim(extracted_shift_vectors_file) > 0) call delete_file(extracted_shift_vectors_file, ierr)
        call delete_file("manifest.txt", ierr)
        
        deallocate(keys, filenames)
    end subroutine read_tox_data

end module tox_archive

!> R-callable generic archive creation with arrays of keys and filenames
subroutine create_zip_archive_generic_R(zip_filename, zip_len, &
                                             keys, keys_len, keys_count, &
                                             filenames, filenames_len, filenames_count, &
                                             ierr)
    use tox_archive, only: create_zip_archive
    use tox_conversions, only: c_char_2d_as_string, c_char_1d_as_string
    use iso_c_binding, only: c_char
    use tox_errors, only: is_err, set_ok, set_err_once, ERR_INVALID_INPUT
    use iso_fortran_env, only: int32
    
    integer(int32), intent(in) :: zip_len
    !! length of the zip file
    character(kind=c_char, len=1), intent(in) :: zip_filename(zip_len)
    !! Zip filename as c_chars
    integer(int32), intent(in) :: keys_count
    !! Number of keys
    integer(int32), intent(in) :: keys_len
    !! Length of the keys
    character(kind=c_char, len=1), intent(in) :: keys(keys_len, keys_count)
    !! Keys as c_chars
    integer(int32), intent(in) :: filenames_count
    !! Number of files  
    integer(int32), intent(in) :: filenames_len
    !! Length of the filenames
    character(kind=c_char, len=1), intent(in) :: filenames(filenames_len, filenames_count)
    !! Filenames as c_chars
    integer(int32), intent(out) :: ierr
    !! Error code
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: f_keys(:)
    character(len=:), allocatable :: f_filenames(:)
    
    call set_ok(ierr)
    
    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if(is_err(ierr)) return
    
    call c_char_2d_as_string(keys, f_keys, ierr)
    if(is_err(ierr)) return
    
    call c_char_2d_as_string(filenames, f_filenames, ierr)
    if(is_err(ierr)) return
    
    ! Validate array sizes
    if (size(f_keys) /= size(f_filenames)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
    
    call create_zip_archive(f_zip_filename, f_keys, f_filenames, ierr)
    
end subroutine create_zip_archive_generic_R

!> Wrapper to extract all files from a zip archive
subroutine extract_zip_archive_generic_R(zip_filename, filename_len, ierr)
    use iso_c_binding, only: c_char, c_null_char
    use iso_fortran_env, only: int32
    use tox_conversions, only: c_char_1d_as_string
    use tox_archive, only: extract_zip_archive
    use tox_errors, only: set_ok, is_err
    implicit none

    integer(int32), intent(in) :: filename_len
    !! Length of the filename
    character(kind=c_char, len=1), intent(in) :: zip_filename(filename_len)
    !! Zip filename as c_chars
    integer(int32), intent(out) :: ierr
    !! Error code
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: keys(:), filenames(:)
    
    call set_ok(ierr)
    
    ! Convert C string to Fortran string
    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if(is_err(ierr)) return
    
    ! Extract the archive - this will extract all files to current directory
    call extract_zip_archive(f_zip_filename, keys, filenames, ierr)
    
    ! Clean up allocated arrays
    if (allocated(keys)) deallocate(keys)
    if (allocated(filenames)) deallocate(filenames)
end subroutine extract_zip_archive_generic_R

!> C binding for generic archive creation with arrays of keys and filenames
subroutine create_zip_archive_generic_c(zip_filename, zip_len, &
                                             keys, keys_len, keys_count, &
                                             filenames, filenames_len, filenames_count, &
                                             ierr) bind(C, name="create_zip_archive_generic_c")
    use tox_archive, only: create_zip_archive
    use tox_conversions, only: c_char_2d_as_string, c_char_1d_as_string
    use iso_c_binding, only: c_int, c_char
    use tox_errors, only: is_err, set_ok, set_err_once, ERR_INVALID_INPUT
    use iso_fortran_env, only: int32
    
    ! Input arguments
    integer(c_int), intent(in), value :: zip_len
    !! Length of the zip filename
    character(kind=c_char, len=1), intent(in) :: zip_filename(zip_len)
    !! Zip Filename as c_chars
    integer(c_int), intent(in), value :: keys_count
    !! number of keys
    integer(c_int), intent(in), value :: keys_len
    !! lengths of the keys
    character(kind=c_char, len=1), intent(in) :: keys(keys_len, keys_count)
    !! Keys as c_chars
    integer(c_int), intent(in), value :: filenames_count
    !! Number of files  
    integer(c_int), intent(in), value :: filenames_len
    !! Length of the filenames
    character(kind=c_char, len=1), intent(in) :: filenames(filenames_len, filenames_count)
    !! Filenames as c_chars
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: f_keys(:)
    character(len=:), allocatable :: f_filenames(:)
    
    call set_ok(ierr)
    
    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if(is_err(ierr)) return
    
    call c_char_2d_as_string(keys, f_keys, ierr)
    if(is_err(ierr)) return
    
    call c_char_2d_as_string(filenames, f_filenames, ierr)
    if(is_err(ierr)) return
    
    ! Validate array sizes
    if (size(f_keys) /= size(f_filenames)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
    
    call create_zip_archive(f_zip_filename, f_keys, f_filenames, ierr)
    
end subroutine create_zip_archive_generic_c

!> C binding for extract_zip_archive - can be called directly from Python via ctypes
subroutine extract_zip_archive_generic_c(zip_filename, filename_len, ierr) &
                                bind(C, name="extract_zip_archive_generic_c")
    use tox_archive, only: extract_zip_archive
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, is_err
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32

    ! Input arguments
    integer(c_int), intent(in), value :: filename_len
    !! Length of the filename
    character(kind=c_char, len=1), intent(in) :: zip_filename(filename_len)
    !! Zip filename length
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: keys(:), filenames(:)
    
    call set_ok(ierr)
    
    ! Convert C string to Fortran string
    call c_char_1d_as_string(zip_filename, f_zip_filename, ierr)
    if(is_err(ierr)) return
    
    call extract_zip_archive(f_zip_filename, keys, filenames, ierr)
    
end subroutine extract_zip_archive_generic_c