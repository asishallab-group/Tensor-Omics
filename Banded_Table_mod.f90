module BandedTable_mod
    use iso_fortran_env, only: real64, character_storage_size, error_unit, iostat_end, iostat_eor
    implicit none

    private ! Default to private, expose specifics

    ! --- Constants ---
    integer, parameter :: TC_INT = 1
    integer, parameter :: TC_REAL = 2
    integer, parameter :: TC_STRING = 3
    character(len=1), parameter :: DEFAULT_SEPARATOR = achar(9) ! Tab (achar is intrinsic)
    integer, parameter :: MAX_LINE_LENGTH = 4096 ! Max expected line length for buffers
    integer, parameter :: MAX_FIELD_LENGTH = 1024 ! Max expected field length for buffers

    ! --- Derived Types ---

    type :: TableColumn
        character(len=:), allocatable :: name
        integer :: column_index = 0  ! Original index in the file
        integer :: type_code = 0     ! 1=int, 2=real, 3=string
        integer :: start_index = 0   ! Starting index within the specific type band (1-based)
                                      ! For strings, this is the index into the offsets array
    end type TableColumn

    type :: Table
        integer :: n_rows = 0
        integer :: n_columns = 0
        ! TARGET attribute removed from here
        type(TableColumn), allocatable :: columns(:)
        ! Data Bands
        integer,          allocatable :: int_values(:)    ! Stores all integer columns contiguously
        real(real64),     allocatable :: real_values(:)   ! Stores all real columns contiguously
        ! Declaration is array of single characters:
        character(len=1), allocatable :: string_chars(:)  ! Stores all string characters contiguously
        integer,          allocatable :: string_offsets(:)! Stores start index (+1 for end) of each string in string_chars
                                                          ! Size: (number of string cells + 1)
    end type Table

    ! --- Public Interface ---
    public :: TableColumn, Table
    public :: parse_table
    public :: print_preview
    public :: get_column_by_name
    public :: destroy_table
    ! Commented out unimplemented procedures:
    ! public :: write_binary, read_binary

contains

    !=======================================================================
    ! destroy_table: Deallocates memory associated with a Table object
    !=======================================================================
    subroutine destroy_table(tbl)
        type(Table), intent(inout) :: tbl
        integer :: j

        ! Compiler should now recognize tbl%columns correctly
        if (allocated(tbl%columns)) then
            do j = 1, size(tbl%columns) ! Use size() on the array component
               if (allocated(tbl%columns(j)%name)) deallocate(tbl%columns(j)%name)
            end do
            deallocate(tbl%columns)
        end if
        if (allocated(tbl%int_values))     deallocate(tbl%int_values)
        if (allocated(tbl%real_values))    deallocate(tbl%real_values)
        if (allocated(tbl%string_chars))   deallocate(tbl%string_chars)
        if (allocated(tbl%string_offsets)) deallocate(tbl%string_offsets)

        tbl%n_rows = 0
        tbl%n_columns = 0
    end subroutine destroy_table

    !=======================================================================
    ! parse_table: Reads a delimited text file into the banded Table structure
    !=======================================================================
    subroutine parse_table(filename, tbl, header, separator, ierror)
        character(len=*), intent(in) :: filename
        type(Table), intent(out)      :: tbl 
        logical, intent(in)           :: header
        character(len=1), intent(in)  :: separator
        integer, intent(out)          :: ierror ! 0 = success, >0 = error code

        ! --- Local variables ---
        integer :: lun, stat, i, j, k, line_count, current_char_pos, char_idx
        integer :: n_int_cols, n_real_cols, n_string_cols, total_string_chars
        integer :: current_int_start, current_real_start, current_string_start_idx
        integer, allocatable :: column_types(:) ! 
        character(len=MAX_LINE_LENGTH) :: line ! Use fixed-size buffer for reading lines
        character(len=:), allocatable :: temp_str_alloc ! For processing fields
        character(len=:), allocatable :: fields(:)
        character(len=MAX_FIELD_LENGTH) :: field_buffer ! Buffer for reading fields safely
        integer :: val_int
        real(real64) :: val_real
        logical :: conversion_ok
        character(len=:), allocatable :: header_names(:) ! Still allocatable, length determined on allocation

        ierror = 0
        call destroy_table(tbl) ! Ensure tbl starts clean

        ! --- File Open ---
        open(newunit=lun, file=trim(filename), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            ierror = 1 ! File open error
            write(error_unit,'(a,a,i0)') 'Error opening file: ', trim(filename), stat
            return
        end if

        ! --- Pass 1: Determine dimensions, column names (if header), and final types ---
        line_count = 0
        total_string_chars = 0
        tbl%n_columns = 0

        ! Read Header (if it exists)
        if (header) then
            read(lun, '(a)', iostat=stat) line
            if (stat /= 0) then
                ierror = 2 ! Error reading header or empty file
                close(lun)
                write(error_unit,'(a)') 'Error reading header line or file is empty.'
                return
            end if
            call split_line(trim(line), separator, fields, stat)
            ! Check stat immediately after split_line
            if (stat /= 0) then
                 ierror = 3 ! Error splitting header (or allocation failed in split_line)
                 close(lun)
                 write(error_unit,'(a)') 'Error splitting header line.'
                 ! fields might not be allocated if split_line failed early
                 if(allocated(fields)) deallocate(fields)
                 return
            end if
            tbl%n_columns = size(fields)
            ! Allocate header_names with specific length
            allocate(character(len=MAX_FIELD_LENGTH) :: header_names(tbl%n_columns), stat=stat)
             if (stat /= 0) then
                ierror = 21 ! Allocation error
                write(error_unit,*) 'Error: Failed to allocate memory for header_names.'
                close(lun)
                if(allocated(fields)) deallocate(fields)
                return
            end if
            do j = 1, tbl%n_columns
                ! Trim before assigning to fixed-length buffer to avoid overflow issues if field is too long
                 if (len_trim(fields(j)) > MAX_FIELD_LENGTH) then
                     write(error_unit,*) 'Warning: Header name truncated: ', trim(fields(j)(1:MAX_FIELD_LENGTH))
                     header_names(j) = trim(fields(j)(1:MAX_FIELD_LENGTH))
                 else
                     header_names(j) = trim(fields(j))
                 end if
            end do
            if(allocated(fields)) deallocate(fields)
        end if

        ! Determine column count (if no header) and initial types from first data line
        allocate(column_types(0)) ! Placeholder
        do
            read(lun, '(a)', iostat=stat) line
            if (stat == iostat_end) exit ! End of file reached before finding data lines
            if (stat /= 0) then
                 ierror = 4 ! Error reading first data line
                 close(lun)
                 write(error_unit,'(a)') 'Error reading first data line.'
                 if(allocated(header_names)) deallocate(header_names)
                 if(allocated(column_types)) deallocate(column_types)
                 return
            end if
            line = trim(line)
            if (len(line) == 0) cycle ! Skip empty lines

            call split_line(line, separator, fields, stat)
            ! Check stat immediately after split_line
            if (stat /= 0) then
                 ierror = 5 ! Error splitting first data line
                 close(lun)
                 write(error_unit,'(a)') 'Error splitting first data line.'
                 if(allocated(header_names)) deallocate(header_names)
                 if(allocated(column_types)) deallocate(column_types)
                 if(allocated(fields)) deallocate(fields)
                 return
            end if

            if (tbl%n_columns == 0) then ! No header, determine n_columns now
                tbl%n_columns = size(fields)
                if (tbl%n_columns == 0) then
                   ierror = 6 ! First data line is empty or has no fields
                   close(lun)
                   write(error_unit,'(a)') 'First data line is empty or has no fields.'
                   if(allocated(header_names)) deallocate(header_names)
                   if(allocated(column_types)) deallocate(column_types)
                   if(allocated(fields)) deallocate(fields)
                   return
                end if
            else ! Header was present, check consistency
                if (size(fields) /= tbl%n_columns) then
                    ierror = 7 ! Column count mismatch
                    close(lun)
                    write(error_unit,'(a,i0,a,i0)') 'Column count mismatch on first data line. Expected ', &
                          tbl%n_columns, ', found ', size(fields)
                    if(allocated(header_names)) deallocate(header_names)
                    if(allocated(column_types)) deallocate(column_types)
                    if(allocated(fields)) deallocate(fields)
                    return
                end if
            end if

            ! Allocate and determine initial types
            if(allocated(column_types)) deallocate(column_types)
            allocate(column_types(tbl%n_columns))
            column_types = 0 ! 0 = undetermined yet
            do j = 1, tbl%n_columns
                temp_str_alloc = trim(fields(j)) ! Use allocatable temp string
                if (len(temp_str_alloc) == 0) then
                   ! Treat empty fields initially as potentially integer, will upgrade if needed
                   column_types(j) = TC_INT
                   cycle
                end if
                ! Try Integer
                read(temp_str_alloc, *, iostat=stat) val_int
                if (stat == 0) then
                    column_types(j) = TC_INT
                    cycle
                end if
                ! Try Real
                read(temp_str_alloc, *, iostat=stat) val_real
                if (stat == 0) then
                    column_types(j) = TC_REAL
                    cycle
                end if
                ! Default to String
                column_types(j) = TC_STRING
                total_string_chars = total_string_chars + len(temp_str_alloc) ! Add length for first row
            end do
            line_count = 1
            if(allocated(fields)) deallocate(fields)
            if(allocated(temp_str_alloc)) deallocate(temp_str_alloc)
            exit ! First data line processed, exit this loop
        end do

        ! Continue reading lines to refine types and count rows/string chars
        do
            read(lun, '(a)', iostat=stat) line
            if (stat == iostat_end) exit ! Normal end of file
            if (stat /= 0) then
                ierror = 8 ! Error reading subsequent data line
                write(error_unit,'(a,i0)') 'Error reading data line ', line_count + 1
                goto 999 ! Cleanup needed
            end if
            line = trim(line)
            if (len(line) == 0) cycle ! Skip empty lines

            call split_line(line, separator, fields, stat)
             ! Check stat immediately after split_line
             if (stat /= 0 .or. size(fields) /= tbl%n_columns) then
                 ierror = 9 ! Error splitting line or column count mismatch
                 write(error_unit,'(a,i0,a,i0)') 'Error splitting line or column count mismatch on line ', &
                       line_count + 1, '. Expected ', tbl%n_columns, ' columns.'
                 if(allocated(fields)) deallocate(fields)
                 goto 999 ! Cleanup needed
            end if

            line_count = line_count + 1
            do j = 1, tbl%n_columns
                temp_str_alloc = trim(fields(j))
                if (len(temp_str_alloc) == 0) cycle ! Empty fields don't change type

                if (column_types(j) == TC_INT) then
                    read(temp_str_alloc, *, iostat=stat) val_int
                    if (stat /= 0) then ! Failed int conversion
                       read(temp_str_alloc, *, iostat=stat) val_real
                       if (stat == 0) then ! Success converting to real
                           column_types(j) = TC_REAL ! Upgrade type
                       else ! Failed real too
                           column_types(j) = TC_STRING ! Upgrade type
                           total_string_chars = total_string_chars + len(temp_str_alloc)
                       end if
                    end if
                else if (column_types(j) == TC_REAL) then
                    read(temp_str_alloc, *, iostat=stat) val_real
                    if (stat /= 0) then ! Failed real conversion
                        column_types(j) = TC_STRING ! Upgrade type
                        total_string_chars = total_string_chars + len(temp_str_alloc)
                    end if
                else if (column_types(j) == TC_STRING) then
                     total_string_chars = total_string_chars + len(temp_str_alloc) !add length
                end if
            end do
            if(allocated(fields)) deallocate(fields)
            if(allocated(temp_str_alloc)) deallocate(temp_str_alloc)
        end do

        tbl%n_rows = line_count
        rewind(lun) ! Go back to the beginning for the second pass

        ! --- Check if any data was read ---
        if (tbl%n_rows == 0 .or. tbl%n_columns == 0) then
            ierror = 10 ! No data rows found
            write(error_unit,'(a)') 'No valid data rows found in the file.'
            goto 999
        end if

        ! --- Allocate Table Structure ---
        n_int_cols = count(column_types == TC_INT)
        n_real_cols = count(column_types == TC_REAL)
        n_string_cols = count(column_types == TC_STRING)

        allocate(tbl%columns(tbl%n_columns), stat=stat)
        if (stat /= 0) then; ierror = 11; goto 998; end if

        if (n_int_cols > 0) then
            allocate(tbl%int_values(n_int_cols * tbl%n_rows), stat=stat)
            if (stat /= 0) then; ierror = 12; goto 998; end if
            tbl%int_values = 0 ! Initialize (optional, good practice)
        end if
        if (n_real_cols > 0) then
            allocate(tbl%real_values(n_real_cols * tbl%n_rows), stat=stat)
             if (stat /= 0) then; ierror = 13; goto 998; end if
            tbl%real_values = 0.0_real64 ! Initialize
        end if
        if (n_string_cols > 0) then
            ! Allocate string_chars as array of len=1 characters
            allocate(tbl%string_chars(max(1, total_string_chars)), stat=stat)
             if (stat /= 0) then; ierror = 14; goto 998; end if
            tbl%string_chars = ' ' ! Initialize
            ! Allocate offsets: number of string cells + 1
            allocate(tbl%string_offsets(n_string_cols * tbl%n_rows + 1), stat=stat)
             if (stat /= 0) then; ierror = 15; goto 998; end if
            tbl%string_offsets = 0 ! Initialize
        end if

        ! Populate column metadata
        current_int_start = 1
        current_real_start = 1
        current_string_start_idx = 1 ! Index within the string *offsets* array
        k = 0 ! Counter for string columns found so far
        do j = 1, tbl%n_columns
            tbl%columns(j)%column_index = j
            tbl%columns(j)%type_code = column_types(j)
            if (header) then
                ! Allocate name with specific length, copy from fixed-length header_names
                allocate(character(len=len_trim(header_names(j))) :: tbl%columns(j)%name, stat=stat)
                 if (stat /= 0) then; ierror = 22; goto 998; end if
                tbl%columns(j)%name = trim(header_names(j))
            else
                ! Increase length for default names
                allocate(character(len=20) :: tbl%columns(j)%name, stat=stat) ! Increased length for "Col_j"
                 if (stat /= 0) then; ierror = 23; goto 998; end if
                ! Internal write should be safer now
                write(tbl%columns(j)%name, '(a,i0)') 'Col_', j
            end if

            select case (tbl%columns(j)%type_code)
            case (TC_INT)
                tbl%columns(j)%start_index = current_int_start
                current_int_start = current_int_start + tbl%n_rows
            case (TC_REAL)
                tbl%columns(j)%start_index = current_real_start
                current_real_start = current_real_start + tbl%n_rows
            case (TC_STRING)
                tbl%columns(j)%start_index = current_string_start_idx
                current_string_start_idx = current_string_start_idx + tbl%n_rows
                k = k + 1
            end select
        end do
        if (allocated(header_names)) deallocate(header_names)
        if (allocated(column_types)) deallocate(column_types)


        ! --- Pass 2: Read data and populate bands ---
        if (header) then
            read(lun, '(a)', iostat=stat) line ! Skip header line
            if (stat /= 0) then; ierror = 16; goto 999; end if ! Should not happen
        end if

        current_char_pos = 0 ! Current end position (0-based index) in string_chars
        if (n_string_cols > 0) tbl%string_offsets(1) = 1 ! First string starts at char 1 (1-based index)

        do i = 1, tbl%n_rows
            read(lun, '(a)', iostat=stat) line
             if (stat == iostat_end .and. i <= tbl%n_rows) then ! Check if EOF happened before expected rows read
                 ierror = 17 ! Premature end of file
                 write(error_unit,*) 'Error: Premature end of file during second pass. Expected ', &
                                     tbl%n_rows, ' data rows, found ', i-1
                 goto 999
             else if (stat /= 0 .and. stat /= iostat_end) then
                 ierror = 18 ! Read error during second pass
                 write(error_unit,*) 'Error: File read error during second pass on data line ', i
                 goto 999
            end if
            line = trim(line)
            if (len(line) == 0) then
                 ! Handle empty line - fill with defaults.
                 do j = 1, tbl%n_columns
                    select case(tbl%columns(j)%type_code)
                    case(TC_INT)
                       tbl%int_values(tbl%columns(j)%start_index + i - 1) = 0
                    case(TC_REAL)
                       tbl%real_values(tbl%columns(j)%start_index + i - 1) = 0.0_real64
                    case(TC_STRING)
                       k = tbl%columns(j)%start_index + i - 1 ! Index into offsets array for this cell
                       ! Empty string: next offset = current offset
                       ! Ensure previous offset was set correctly (should be okay if k > 0)
                       if (k == 0) then ! First string cell overall needs offset(1) initialized
                          tbl%string_offsets(1) = 1
                       end if
                       tbl%string_offsets(k + 1) = tbl%string_offsets(k)
                    end select
                 end do
                 cycle
            end if

            call split_line(line, separator, fields, stat)
            ! Check stat immediately after split_line
            if (stat /= 0 .or. size(fields) /= tbl%n_columns) then
                ierror = 19 ! Split error or column mismatch on pass 2
                write(error_unit,*) 'Error: Splitting error or column mismatch on pass 2, line ', i
                if(allocated(fields)) deallocate(fields)
                goto 999
            end if

            do j = 1, tbl%n_columns
                ! Use fixed-size buffer for read safety, but process trimmed version
                field_buffer = fields(j)
                temp_str_alloc = trim(field_buffer)

                select case (tbl%columns(j)%type_code)
                case (TC_INT)
                    if (len(temp_str_alloc) > 0) then
                       read(temp_str_alloc, *, iostat=stat) val_int
                       if (stat /= 0) then
                           write(error_unit,*) 'Warning: Int conversion failed on pass 2 for [', temp_str_alloc, &
                                               '] at row ', i, ', col ', j, '. Storing 0.'
                           val_int = 0
                       end if
                    else
                       val_int = 0 ! Default for empty int field
                    end if
                    tbl%int_values(tbl%columns(j)%start_index + i - 1) = val_int

                case (TC_REAL)
                     if (len(temp_str_alloc) > 0) then
                        read(temp_str_alloc, *, iostat=stat) val_real
                        if (stat /= 0) then
                             write(error_unit,*) 'Warning: Real conversion failed on pass 2 for [', temp_str_alloc, &
                                                '] at row ', i, ', col ', j, '. Storing 0.0.'
                            val_real = 0.0_real64
                        end if
                     else
                        val_real = 0.0_real64 ! Default for empty real field
                     end if
                    tbl%real_values(tbl%columns(j)%start_index + i - 1) = val_real

                case (TC_STRING)
                    k = tbl%columns(j)%start_index + i - 1 ! Index into offsets array for this cell
                    ! Current string data starts at index tbl%string_offsets(k)
                    ! Set the *next* offset first based on the length of the current string
                    tbl%string_offsets(k + 1) = tbl%string_offsets(k) + len(temp_str_alloc)

                    ! Now copy the characters if the string is not empty
                    if (len(temp_str_alloc) > 0) then
                       ! Check bounds before copying
                       if (tbl%string_offsets(k+1) - 1 > size(tbl%string_chars)) then
                           ierror = 20 ! String buffer overflow
                           write(error_unit,*) 'Error: String character buffer overflow detected during copy.'
                           if(allocated(fields)) deallocate(fields)
                           if(allocated(temp_str_alloc)) deallocate(temp_str_alloc)
                           goto 999
                       end if
                       ! *** Correction: Copy characters explicitly ***
                       do char_idx = 1, len(temp_str_alloc)
                           tbl%string_chars(tbl%string_offsets(k) + char_idx - 1) = temp_str_alloc(char_idx:char_idx)
                       end do
                    end if
                    ! Update overall character position tracker (optional, for final check)
                    current_char_pos = tbl%string_offsets(k+1) - 1

                end select
            end do ! End loop columns (j)
            if(allocated(fields)) deallocate(fields)
            if(allocated(temp_str_alloc)) deallocate(temp_str_alloc)
        end do ! End loop over rows (i)

        ! --- Final checks ---
        if (n_string_cols > 0 .and. current_char_pos /= total_string_chars) then
             write(error_unit,*) 'Warning: Final character count position (', current_char_pos, &
                                 ') differs from Pass 1 estimate (', total_string_chars, '). Check data/parsing.'
             ! This might indicate issues with line endings or special chars counted differently.
        end if

999 continue ! Target for cleanup after error during passes
        close(lun)
        ! Deallocate temporary arrays even if error occurred during loops
        if (allocated(header_names)) deallocate(header_names)
        if (allocated(column_types)) deallocate(column_types)
        if (allocated(fields)) deallocate(fields) ! Ensure fields is deallocated if error happened mid-loop
        if (allocated(temp_str_alloc)) deallocate(temp_str_alloc)
        ! If error occurred, ensure the table is deallocated
        if (ierror /= 0) then
            call destroy_table(tbl)
        end if
        return

998 continue ! Target for cleanup after allocation error
        write(error_unit,'(a,i0)') 'Error: Failed to allocate memory. Error code: ', ierror
        close(lun, iostat=stat) ! Close file even if allocation failed
        if (allocated(header_names)) deallocate(header_names)
        if (allocated(column_types)) deallocate(column_types)
        call destroy_table(tbl) ! Clean up partially allocated table
        return
    end subroutine parse_table

    !=======================================================================
    ! split_line: Helper to split a string by a separator
    ! Allocates the output array 'parts'. Caller must deallocate.
    !=======================================================================
    subroutine split_line(line, separator, parts, stat)
        character(len=*), intent(in) :: line
        character(len=1), intent(in) :: separator
        character(len=:), allocatable, intent(out) :: parts(:)
        integer, intent(out) :: stat ! 0 = ok, 1 = allocation error

        character(len=len(line)) :: work_line
        integer :: n_parts, i, start_pos, current_alloc
        character(len=:), allocatable :: temp_parts(:)
        integer :: count
        integer :: max_part_len ! To determine length for reallocation

        stat = 0
        work_line = trim(line)
        count = 0

        if (len(work_line) == 0) then
           if (allocated(parts)) deallocate(parts)
           ! Ensure output stat is set on allocation failure
           allocate(character(len=0) :: parts(0), stat=stat) ! Use the output stat
           if (stat /= 0) then
               write(error_unit,*) 'Split_line: Error allocating zero-sized parts array.'
               ! stat is already non-zero from allocate
           end if
           return ! Return, caller will check stat
        end if

        ! First pass: count parts to allocate correctly
        n_parts = 1
        do i = 1, len(work_line)
           if (work_line(i:i) == separator) then
              n_parts = n_parts + 1
           end if
        end do

        if (allocated(parts)) deallocate(parts)
        ! Allocate parts with explicit length (max possible needed)
        allocate(character(len=len(line)) :: parts(n_parts), stat=stat)
        if (stat /= 0) then
            write(error_unit,*) 'Split_line: Error allocating parts array.'
            return ! Allocation error (stat is already non-zero)
        end if

        ! Second pass: extract parts
        count = 0
        start_pos = 1
        do i = 1, len(work_line)
           if (work_line(i:i) == separator) then
              count = count + 1
              parts(count) = work_line(start_pos:i-1)
              start_pos = i + 1
           end if
        end do
        ! Add the last part
        count = count + 1
        parts(count) = work_line(start_pos:)

        ! Trim lengths (optional but good practice)
        ! Allocate temp_parts with explicit length before move_alloc
        allocate(character(len=len(line)) :: temp_parts(n_parts), stat=stat)
         if (stat /= 0) then
            write(error_unit,*) 'Split_line: Error allocating temp_parts array.'
            ! Cannot proceed with trimming, leave 'parts' as is (with max length)
            ! stat is already non-zero
            return
        end if
        max_part_len = 0
        do i = 1, n_parts
            temp_parts(i) = trim(parts(i))
            max_part_len = max(max_part_len, len(temp_parts(i)))
        end do
        ! Now reallocate parts to the actual needed max length
        if (allocated(parts)) deallocate(parts)
        allocate(character(len=max_part_len) :: parts(n_parts), stat=stat)
         if (stat /= 0) then
            write(error_unit,*) 'Split_line: Error re-allocating parts array with trimmed length.'
             if (allocated(temp_parts)) deallocate(temp_parts)
             ! stat is already non-zero
            return
        end if
        parts = temp_parts ! Copy trimmed strings
        if (allocated(temp_parts)) deallocate(temp_parts)

    end subroutine split_line


    !=======================================================================
    ! get_column_by_name: Find column metadata by name (case-sensitive)
    !=======================================================================
    function get_column_by_name(tbl, name) result(col)
        ! Add TARGET attribute to the dummy argument tbl
        type(Table), intent(in), target :: tbl
        character(len=*), intent(in) :: name
        type(TableColumn), pointer :: col

        integer :: j
        logical :: found

        col => null() ! Return null pointer if not found
        if (.not. allocated(tbl%columns)) return

        found = .false.
        do j = 1, tbl%n_columns
            ! Pointer assignment should now work because tbl is a target
            if (trim(tbl%columns(j)%name) == trim(name)) then
                col => tbl%columns(j)
                found = .true.
                exit
            end if
        end do

    end function get_column_by_name

    !=======================================================================
    ! print_preview: Prints the first few rows of the table
    !=======================================================================
    subroutine print_preview(tbl, max_rows, max_cols)
        ! Add TARGET attribute to the dummy argument tbl
        type(Table), intent(in), target :: tbl
        integer, intent(in), optional :: max_rows
        integer, intent(in), optional :: max_cols

        integer :: r_max, c_max, i, j, k, start_char, end_char, len_str, char_idx
        character(len=25) :: fmt_header, fmt_int, fmt_real, fmt_str, buffer
        character(len=:), allocatable :: temp_string ! Use allocatable for string extraction
        type(TableColumn), pointer :: p_col
        integer :: stat ! Declared locally

        if (tbl%n_rows == 0 .or. tbl%n_columns == 0) then
            print *, "Table is empty."
            return
        end if

        r_max = 5
        if (present(max_rows)) r_max = max_rows
        r_max = min(r_max, tbl%n_rows)

        c_max = 10
        if (present(max_cols)) c_max = max_cols
        c_max = min(c_max, tbl%n_columns)

        ! Print header
        write(fmt_header, '(a,i0,a)') '(a15,', c_max, '(1x,a15))'
        write(*, fmt_header) 'Row', (trim(adjustl(tbl%columns(j)%name)), j=1, c_max) ! adjustl for alignment
        write(*,'(a)') repeat('-', 15 + c_max * 16)

        ! Print data rows
        do i = 1, r_max
            ! Use write with advance='no'
            write(*,'(i15)', advance='no') i ! Row number

            do j = 1, c_max
                ! Pointer assignment should now work because tbl is a target
                p_col => tbl%columns(j)
                select case (p_col%type_code)
                case(TC_INT)
                    write(buffer,'(i15)') tbl%int_values(p_col%start_index + i - 1)
                case(TC_REAL)
                    write(buffer,'(es15.7)') tbl%real_values(p_col%start_index + i - 1)
                case(TC_STRING)
                    k = p_col%start_index + i - 1 ! Index in offsets array
                    start_char = tbl%string_offsets(k)
                    end_char   = tbl%string_offsets(k+1) - 1
                    len_str = end_char - start_char + 1
                    if (len_str > 0) then
                       ! Explicitly construct scalar string from char array slice
                       allocate(character(len=len_str) :: temp_string, stat=stat)
                       if (stat /= 0) then
                           buffer = '*ALLOC ERR*'
                       else
                           ! Copy characters one by one from the array slice
                           do char_idx = 1, len_str
                               temp_string(char_idx:char_idx) = tbl%string_chars(start_char + char_idx - 1)
                           end do
                           ! Assign constructed string to fixed buffer (truncates/pads)
                           buffer = temp_string
                           deallocate(temp_string)
                       end if
                    else
                       buffer = '' ! Empty string
                    end if
                    ! Ensure buffer is formatted correctly for printing (left-aligned)
                    buffer = adjustl(buffer)

                case default
                    buffer = '?' ! Unknown type
                    buffer = adjustl(buffer)
                end select
                 ! Use write with advance='no'
                 write(*, '(1x, a15)', advance='no') trim(buffer)
            end do ! End loop columns (j)
            write(*,*) ! Newline after each row
        end do ! End loop rows (i)

        if (r_max < tbl%n_rows .or. c_max < tbl%n_columns) then
             print *, '... (preview truncated) ...'
             print *, 'Total rows:', tbl%n_rows, ', Total columns:', tbl%n_columns
        end if
        ! Deallocate temp_string if it was allocated in the loop and loop exited early (though unlikely here)
        if (allocated(temp_string)) deallocate(temp_string)

    end subroutine print_preview

end module BandedTable_mod
