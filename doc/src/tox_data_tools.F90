#include "../macros.h"

module tox_data_tools
    use safeguard
    use iso_fortran_env, only: real64, int32, iostat_end
    use tox_errors, only: set_ok, set_err_once, is_err, check_io_stat
    use tox_errors, only: ERR_INVALID_INPUT, ERR_FILE_OPEN, ERR_READ_DATA
    use f42_array_utils, only :check_okay_ioerror, ERR_SIZE_MISMATCH
    use config, only: DEBUG
    implicit none
    private

    public :: read_gene_ids_from_tsv_file
    public :: read_expression_vectors_tsv
    public :: read_orthofinder_file
    public :: split_string
    public :: get_unassigned_mask, apply_unassigned_mask

    character(len=*), parameter :: DELIMS = ', ' // char(9)

contains

!> Read expression vectors from csv/tsv files
subroutine read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, &
                             n_header_rows, gene_col, value_cols, start_row, ierr, delimiter)
    use f42_xxh3_hashmap
    use ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: file_list(:)
    !! List of files to read from
    character(len=*), intent(in) :: gene_ids(:)
    !! Array of gene IDS
    real(real64), intent(inout) :: expression_vectors(:,:)
    !! Array of expression vectors
    integer(int32), intent(in) :: n_header_rows
    !! Number of header rows to skip
    integer(int32), intent(in) :: gene_col
    !! Index of column with gene_ids
    integer(int32), intent(in) :: value_cols(:)
    !! Indicies of columns containing values
    integer(int32), intent(in) :: start_row
    !! Row in the expression vectors to start in
    integer(int32), intent(out) :: ierr
    !! Error code
    character(len=1), intent(in), optional :: delimiter
    !! optional delimiter, default is tab

    integer(int32) :: i, j, k, unit, ios, idx, n_genes, expected_idx, n_value_cols
    integer(int32) :: current_sample, n_columns_in_file, n_valid_cols
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:), test_fields(:)
    integer(int32), allocatable :: valid_cols(:)
    real(real64) :: value
    character(len=len(gene_ids)) :: gene
    character(len=1) :: actual_delimiter
    integer(int32) :: current_row

    ! Hashmap for gene lookup
    type(hashmap_type) :: gene_map

    call set_ok(ierr)
    call set_ok(ios)
    n_genes = size(gene_ids)
    n_value_cols = size(value_cols)
    
    ! Set the delimiter (use tab as default if not provided)
    if (present(delimiter)) then
        actual_delimiter = delimiter
    else
        actual_delimiter = char(9)
    end if
    
    if (n_value_cols <= 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        if(DEBUG) write(*,*) 'Error: No value columns specified.'
        return
    end if

    if (size(gene_ids) < n_genes) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        if(DEBUG) write(*,*) 'Error: gene_ids array is too small.'
        return
    end if

    ! Allocate temporary array for valid columns
    allocate(valid_cols(n_value_cols), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    ! Build gene hashmap for fast lookup
    call hashmap_create(gene_map, initial_size=n_genes)
    do i = 1, n_genes
        call hashmap_put(gene_map, trim(gene_ids(i)), i)
    end do
    
    current_sample = start_row - 1
    
    do i = 1, size(file_list)
        if(DEBUG) write(*,*) 'Reading file: ', trim(file_list(i))
        open(newunit=unit, file=trim(file_list(i)), status='old', action='read', iostat=ios)
        call check_okay_ioerror(ios, ierr, ERR_FILE_OPEN, unit)
        if(is_err(ierr)) then
            call hashmap_destroy(gene_map)
            return
        end if

        ! Skip header rows
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if
        end do

        ! Read first data line to determine number of columns in this file
        read(unit, '(A)', iostat=ios) line
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if (is_err(ierr)) then
            call hashmap_destroy(gene_map)
            return
        end if
        
        call split_string(line, test_fields, ierr, actual_delimiter)
        n_columns_in_file = size(test_fields)
        
        ! Rewind to beginning of data section
        rewind(unit)
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if
        end do

        ! Determine valid columns for this file
        n_valid_cols = 0
        do k = 1, n_value_cols
            if (value_cols(k) <= n_columns_in_file) then
                n_valid_cols = n_valid_cols + 1
                valid_cols(n_valid_cols) = value_cols(k)
            else 
                call set_err_once(ierr, ERR_INVALID_INPUT)
                return 
            end if
        end do
        if (n_valid_cols == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            write(*,*) 'Error: No valid value columns found in file: ', trim(file_list(i))
            close(unit)
            cycle
        end if

        current_row = 0
        do
            if(current_row > size(gene_ids)) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if(DEBUG) write(*,*) 'Provided file contains more lines then expected'
                close(unit)
                return
            end if
            current_row = current_row + 1
            read(unit, '(A)', iostat=ios) line
            if(ios == iostat_end) exit !End of file
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if

            call split_string(line, fields, ierr, actual_delimiter)
            if (size(fields) < max(gene_col, maxval(valid_cols(1:n_valid_cols)))) then
                if(DEBUG) print *, 'Invalid input:'
                if(DEBUG) print *, 'size: ', size(fields)
                if(DEBUG) print *, 'row: ', current_row
                call set_err_once(ierr, ERR_INVALID_INPUT)
                call hashmap_destroy(gene_map)
                return 
            end if

            gene = trim(adjustl(fields(gene_col)))

            ! Use hashmap for gene lookup
            idx = hashmap_get(gene_map, gene)
            if (idx == -1) then
                write(*,*) 'Warning: Gene ', trim(gene), ' not found in master gene list'
                cycle
            end if
            
            ! Read all valid value columns for this gene 
            do k = 1, n_valid_cols
                read(fields(valid_cols(k)), *, iostat=ios) value
                if (ios == 0) then
                    if(ieee_is_finite(value)) then
                        expression_vectors(current_sample + k, idx) = value
                    else
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Non-finite value encountered at row ', current_row, ' column ', valid_cols(k)
                        call hashmap_destroy(gene_map)
                        close(unit)
                        return
                    end if
                else
                    call set_err_once(ierr, ERR_READ_DATA)
                    call hashmap_destroy(gene_map)
                    close(unit)
                    return
                end if
            end do
        end do
        close(unit)
        
        current_sample = current_sample + n_valid_cols

    end do

    call hashmap_destroy(gene_map)
    deallocate(valid_cols)
end subroutine read_expression_vectors_tsv

!> Only read the gene ids from a tsv file
subroutine read_gene_ids_from_tsv_file(filename, gene_ids, n_header_rows, gene_col, ierr)
    character(len=*), intent(in) :: filename
        !! Name of the file
    character(len=*), intent(out) :: gene_ids(:)
        !! gene ids array
    integer(int32), intent(in) :: n_header_rows
        !! number of headers to skip
    integer(int32), intent(in) :: gene_col
        !! Index of the column containing gene ids
    integer(int32), intent(out) :: ierr
        !! Error code

    integer(int32) :: unit, ios, j, row_count
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:)

    call set_ok(ierr)
    call set_ok(ios)
    row_count = 0

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    call check_okay_ioerror(ios, ierr, ERR_FILE_OPEN, unit)
    if(is_err(ierr)) return

    ! Skip header rows
    do j = 1, n_header_rows
        read(unit, '(A)', iostat=ios) line
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) return
    end do

    ! Read data rows
    do
        read(unit, '(A)', iostat=ios) line
        if (ios == iostat_end) exit
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) return
        
        row_count = row_count + 1
        if (row_count > size(gene_ids)) then
            if(DEBUG) print *, 'Row: ', row_count
            if(DEBUG) print *, 'Size: ', size(gene_ids)
            call set_err_once(ierr, ERR_SIZE_MISMATCH)  ! More genes than allocated space
            close(unit)
            return
        end if

        call split_string(line, fields, ierr, char(9))
        if (size(fields) < gene_col) cycle

        gene_ids(row_count) = trim(adjustl(fields(gene_col)))
    end do

    close(unit)
end subroutine read_gene_ids_from_tsv_file

!> Read a family file (Orthofinder)
subroutine read_orthofinder_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    use f42_xxh3_hashmap
    character(len=*), intent(in) :: filename
        !! Name of the file
    character(len=*), intent(in) :: gene_ids(:)
        !! gene ids array
    character(len=*), intent(out) :: family_ids(:)
        !! family ids array
    integer(int32), intent(out) :: gene_to_fam(:)
        !! gene to family mapping
    integer(int32), intent(out) :: ierr
        !! Error code

    integer(int32) :: unit, ios, i, j, fam_idx, gene_idx, n_families, n_genes
    integer(int32) :: pos, start_pos, end_pos 
    character(len=4096) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    character(len=len(family_ids)) :: current_family
    type(hashmap_type) :: gene_map
    integer(int32) :: hashmap_size

    call set_ok(ierr)
    gene_to_fam = 0
    n_families = size(family_ids)
    n_genes = size(gene_ids)
    
    ! Hashmap for efficient lookups
    call hashmap_create(gene_map, initial_size = n_genes)
    
    ! Fill hashmap
    do i = 1, n_genes
        call hashmap_put(gene_map, gene_ids(i), i)
    end do

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (is_err(ios)) then
        if(DEBUG) write(*,*) 'Error opening file: ', trim(filename)
        call set_err_once(ierr, ERR_FILE_OPEN)
        call hashmap_destroy(gene_map)
        return
    end if

    ! skip header
    read(unit, '(A)', iostat=ios) line
    call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
    if(is_err(ierr)) then
        call hashmap_destroy(gene_map)
        RETURN
    end if
    
    fam_idx = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios == iostat_end) exit
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) then
            call hashmap_destroy(gene_map)
            RETURN
        end if

        fam_idx = fam_idx + 1

        call split_string(line, fields, ierr, char(9))


        if (size(fields) < 2) then
            if(DEBUG) write(*,*) 'Warning: Skipping invalid line (less than 2 fields).'
            cycle
        end if

        current_family = trim(adjustl(fields(1)))

        family_ids(fam_idx) = current_family

        ! Process all gene columns with Hashmap-Lookup
        do i = 2, size(fields)
            call split_string(fields(i), genes, ierr, ',')
            if(is_err(ierr)) return
            do j = 1, size(genes)
                if (len_trim(genes(j)) == 0) cycle
                
                ! Hashmap-Lookup
                gene_idx = hashmap_get(gene_map, trim(adjustl(genes(j))))
                if (gene_idx > 0) then
                    gene_to_fam(gene_idx) = fam_idx
                else
                    if(DEBUG) write(*,*) 'Gene not found in hashmap: ', trim(adjustl(genes(j)))
                end if
            end do
        end do
    end do
    close(unit)
    
    call hashmap_destroy(gene_map)
end subroutine read_orthofinder_file

!> Helper subroutine to split strings
subroutine split_string(input, output, ierr, delimiter)
    character(len=*), intent(in) :: input
        !! Input string to split
    character(len=:), allocatable, intent(out) :: output(:)
        !! Output of parts
    character(len=1), optional, intent(in) :: delimiter
        !! delimiter, default is space
    integer(int32), intent(out) :: ierr
        !! Error code
    integer(int32) :: ios
    
    character(len=1) :: delim
    integer(int32) :: n, i, start_pos, end_pos, field_count
    integer(int32), allocatable :: field_starts(:), field_ends(:)
    
    call set_ok(ierr)
    call set_ok(ios)

    if (present(delimiter)) then
        delim = delimiter
    else
        delim = ' '
    end if

    ! First pass: find field boundaries
    allocate(field_starts(len_trim(input)), stat=ios)
    call check_io_stat(ios, ierr)
    allocate(field_ends(len_trim(input)), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return
    
    n = 0
    start_pos = 1
    do i = 1, len_trim(input)
        if (input(i:i) == delim) then
            if (start_pos <= i-1) then
                n = n + 1
                field_starts(n) = start_pos
                field_ends(n) = i-1
            end if
            start_pos = i + 1
        end if
    end do
    
    ! Add last field if exists
    if (start_pos <= len_trim(input)) then
        n = n + 1
        field_starts(n) = start_pos
        field_ends(n) = len_trim(input)
    end if

    ! Allocate output array
    allocate(character(len=maxval(field_ends(1:n) - field_starts(1:n) + 1)) :: output(n), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    ! Extract fields
    do i = 1, n
        output(i) = trim(adjustl(input(field_starts(i):field_ends(i))))
    end do
    
end subroutine split_string

!> Helper to create a mask of genes that are unassigned
subroutine get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
    use iso_fortran_env, only: int32
    implicit none
    
    integer(int32), intent(in) :: gene_to_fam(:)
        !! gene to family mapping
    logical, intent(out) :: mask(size(gene_to_fam))
        !! mask for mapping
    integer(int32), intent(out) :: n_genes_kept
        !! number of genes kept
    
    integer(int32) :: i
    
    n_genes_kept = 0
    do i = 1, size(gene_to_fam)
        mask(i) = (gene_to_fam(i) >= 1)
        if (mask(i)) n_genes_kept = n_genes_kept + 1
    end do
end subroutine get_unassigned_mask

!> Applies a precomputed mask to filter out unassigned genes
subroutine apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)

    character(len=*), allocatable, intent(inout) :: gene_ids(:)
    !! gene ids array
    real(real64), allocatable, intent(inout) :: expression_vectors(:,:)
    !! expression vectors array
    integer(int32), allocatable, intent(inout) :: gene_to_fam(:)
    !! gene to family mapping
    logical, intent(in) :: mask(:)
    !! mask for mapping
    integer(int32), intent(out) :: n_genes_kept
    !! number of genes kept
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: i, j, n_genes_total, n_samples, ios
    character(len=len(gene_ids)), allocatable :: temp_gene_ids(:)
    real(real64), allocatable :: temp_expression_vectors(:,:)
    integer(int32), allocatable :: temp_gene_to_fam(:)

    call set_ok(ierr)
    n_genes_total = size(gene_to_fam)
    n_samples = size(expression_vectors, 1)
    n_genes_kept = count(mask)

    if (n_genes_kept == 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    allocate(temp_gene_ids(n_genes_kept), stat=ios)
    allocate(temp_expression_vectors(n_samples, n_genes_kept), stat=ios)
    allocate(temp_gene_to_fam(n_genes_kept), stat=ios)
    if (ios /= 0) then
        ierr = ios
        return
    end if

    j = 1
    do i = 1, n_genes_total
        if (mask(i)) then
            temp_gene_ids(j) = gene_ids(i)
            temp_expression_vectors(:, j) = expression_vectors(:, i)
            temp_gene_to_fam(j) = gene_to_fam(i)
            j = j + 1
        end if
    end do

    deallocate(gene_ids, expression_vectors, gene_to_fam)
    call move_alloc(temp_gene_ids, gene_ids)
    call move_alloc(temp_expression_vectors, expression_vectors)
    call move_alloc(temp_gene_to_fam, gene_to_fam)
end subroutine


end module tox_data_tools



!> C binding for reading gene IDs from a gene expression tsv file
subroutine read_gene_ids_from_tsv_file_C(filename_raw, fn_len, gene_ids_raw, gene_ids_len, n_genes, &
                                 n_header_rows, gene_col, ierr) bind(C, name="read_gene_ids_from_tsv_file_C")
    use iso_c_binding, only: c_int, c_char
    use tox_errors, only: set_ok, is_err
    use tox_conversions, only: c_char_1d_as_string, string_as_c_char_1d
    use tox_data_tools, only: read_gene_ids_from_tsv_file
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: fn_len       
        !! Length of filename
    character(kind=c_char, len=1), intent(in), target :: filename_raw(fn_len)
        !! Pointer to filename array
    integer(c_int), intent(in), target :: gene_ids_len 
        !! Length of each gene ID string
    integer(c_int), intent(in), target :: n_genes      
        !! Number of genes    
    character(kind=c_char, len=1), intent(out), target :: gene_ids_raw(gene_ids_len, n_genes)
        !! Pointer to gene_ids array
    integer(c_int), intent(in), target :: n_header_rows
        !! Number of header rows to skip
    integer(c_int), intent(in), target :: gene_col
        !! Index of the gene column
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    integer(c_int) :: i

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(gene_ids_len)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_header_rows)
    M_CHECK_NON_NULL(gene_col)
    M_CHECK_NON_NULL(filename_raw)
    M_CHECK_NON_NULL(gene_ids_raw)

    call set_ok(ierr)

    ! Convert filename from raw bytes
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if(is_err(ierr)) return

    call read_gene_ids_from_tsv_file(filename, gene_ids, n_header_rows, gene_col, ierr)
    if(is_err(ierr)) return

    ! Convert gene IDs back to raw bytes for output
    do i = 1, n_genes
        call string_as_c_char_1d(trim(gene_ids(i)), gene_ids_raw(:, i))
    end do
end subroutine read_gene_ids_from_tsv_file_C

!> C binding for reading expression vectors from files
subroutine read_expression_vectors_tsv_C(file_list_raw, file_list_len, n_files, &
                                 gene_ids_raw, gene_ids_len, n_genes, &
                                 expression_vectors, n_samples, &
                                 n_header_rows, gene_col, value_cols, &
                                 n_value_cols, ierr, delimiter_raw) bind(C, name="read_expression_vectors_tsv_C")
    use iso_c_binding, only: c_int, c_double, c_char
    use tox_data_tools, only: read_expression_vectors_tsv
    use tox_errors, only: set_ok, is_err, check_io_stat
    use tox_conversions, only: c_char_1d_as_string, string_as_c_char_1d
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: file_list_len     
        !! Length of each filename
    integer(c_int), intent(in), target :: n_files           
        !! Number of files
    character(kind=c_char, len=1), intent(in), target :: file_list_raw(file_list_len, n_files)    
        !! Pointer to file_list array
    integer(c_int), intent(in), target :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), intent(in), target :: n_genes           
        !! Number of genes
    character(kind=c_char, len=1), intent(in), target :: gene_ids_raw(gene_ids_len, n_genes)     
        !! Pointer to gene_ids array
    integer(c_int), intent(in), target :: n_samples
        !! Number of samples    
    real(c_double), intent(out), target :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors (flat array)
    integer(c_int), intent(in), target :: n_header_rows
        !! Number of header rows
    integer(c_int), intent(in), target :: gene_col
        !! Index of column containing gene ids
    integer(c_int), intent(in), target :: n_value_cols
        !! Number of cols containing values    
    integer(c_int), intent(in), target :: value_cols (n_value_cols)          
        !! Pointer to value_cols array
    integer(c_int), intent(out), target :: ierr
        !! Error code
    character(kind=c_char, len=1), intent(in), target :: delimiter_raw(1)
        !! Delimiter
    integer(c_int) :: start_row
    
    character(len=file_list_len), allocatable :: file_list(:)
    character(len=gene_ids_len), allocatable :: gene_ids(:)
    character(len=:), allocatable :: delimiter
    character(len=:), allocatable :: tmp_str
    integer(c_int) :: i, j, ios

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(file_list_len)
    M_CHECK_NON_NULL(n_files)
    M_CHECK_NON_NULL(gene_ids_len)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_header_rows)
    M_CHECK_NON_NULL(gene_col)
    M_CHECK_NON_NULL(n_value_cols)
    M_CHECK_NON_NULL(file_list_raw)
    M_CHECK_NON_NULL(gene_ids_raw)
    M_CHECK_NON_NULL(expression_vectors)
    M_CHECK_NON_NULL(value_cols)
    M_CHECK_NON_NULL(delimiter_raw)

    start_row = 1
    call set_ok(ierr)

    allocate(file_list(n_files), stat=ios)
    call check_io_stat(ios, ierr)
    allocate(gene_ids(n_genes), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return
   
    ! Convert 2D raw arrays to string arrays
    do i = 1, n_files
      call c_char_1d_as_string(file_list_raw(:, i), tmp_str, ierr)
      if(is_err(ierr)) return
      file_list(i) = trim(tmp_str)
    end do
    
    do i = 1, n_genes
      call c_char_1d_as_string(gene_ids_raw(:, i), tmp_str, ierr)
      if(is_err(ierr)) return
      gene_ids(i) = trim(tmp_str)
    end do

    call read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, &
                                n_header_rows, gene_col, value_cols, &
                                start_row, ierr, delimiter)

    if(is_err(ierr)) return

end subroutine read_expression_vectors_tsv_C

!> C binding for reading family file
subroutine read_orthofinder_file_C(filename_raw, fn_len, gene_ids_raw, gene_ids_len, n_genes, &
                             family_ids_raw, family_ids_len, n_families, gene_to_fam, ierr) bind(C, name="read_orthofinder_file_C")
    use iso_c_binding, only: c_int, c_char
    use tox_errors, only: set_ok, is_err
    use tox_data_tools, only: read_orthofinder_file
    use tox_conversions, only: c_char_1d_as_string, string_as_c_char_1d
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: fn_len            
        !! Length of filename
    character(kind=c_char, len=1), intent(in), target :: filename_raw(fn_len)       
        !! Pointer to filename array
    integer(c_int), intent(in), target :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), intent(in), target :: n_genes           
        !! Number of genes    
    character(kind=c_char, len=1), intent(in), target :: gene_ids_raw(gene_ids_len, n_genes)       
        !! Pointer to gene_ids array
    integer(c_int), intent(in), target :: family_ids_len    
        !! Length of each family ID
    integer(c_int), intent(in), target :: n_families        
        !! Number of families
    character(kind=c_char, len=1), intent(out), target :: family_ids_raw(family_ids_len, n_families)     
        !! Pointer to family_ids array
    integer(c_int), intent(out), target :: gene_to_fam(n_genes)          
        !! Pointer to gene_to_fam array
    integer(c_int), intent(out), target :: ierr
        !! error code
    
    character(len=:), allocatable :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    character(len=family_ids_len) :: family_ids(n_families)
    character(len=:), allocatable :: temp_str
    integer(c_int) :: i

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(gene_ids_len)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(family_ids_len)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(filename_raw)
    M_CHECK_NON_NULL(gene_ids_raw)
    M_CHECK_NON_NULL(family_ids_raw)
    M_CHECK_NON_NULL(gene_to_fam)

    call set_ok(ierr)

    ! Initialize family_ids with spaces
    family_ids = ' '
    
    ! Convert filename from raw bytes
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if(is_err(ierr)) return
    
    ! Convert gene IDs from raw bytes to strings
    do i = 1, n_genes
        call c_char_1d_as_string(gene_ids_raw(:, i), temp_str, ierr)
        if(is_err(ierr)) return
        gene_ids(i) = trim(adjustl(temp_str))
    end do
    
    call read_orthofinder_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    if(is_err(ierr)) return
    
    ! Convert family IDs to raw bytes
    do i = 1, n_families
        call string_as_c_char_1d(trim(adjustl(family_ids(i))), family_ids_raw(:, i))
    end do
end subroutine read_orthofinder_file_C

!> C binding for filtering unassigned genes
subroutine filter_unassigned_genes_C(gene_ids_raw, gene_ids_len, n_genes, &
                                    gene_to_fam, mask, n_genes_kept, ierr) bind(C, name="filter_unassigned_genes_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT, is_err, check_io_stat
    use tox_conversions, only: c_int_as_logical, logical_as_c_int
    use tox_data_tools, only: get_unassigned_mask
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), intent(in), target :: n_genes           
        !! Number of genes
    character(kind=c_char, len=1), intent(in), target :: gene_ids_raw(gene_ids_len, n_genes)       
        !! Pointer to gene_ids array
    integer(c_int), intent(in), target :: gene_to_fam(n_genes)          
        !! Pointer to gene_to_fam array
    integer(c_int), intent(out), target :: mask(n_genes)                 
        !! Pointer to mask array
    integer(c_int), intent(out) :: n_genes_kept
        !! number of genes kept
    integer(c_int), intent(out), target :: ierr
        !! Error code
    
    integer(int32) :: i, ios
    logical, allocatable :: mask_logical(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(gene_ids_len)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(gene_ids_raw)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(mask)

    call set_ok(ierr)
    call set_ok(ios)

    allocate(mask_logical(n_genes), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    call get_unassigned_mask(gene_to_fam, mask_logical, n_genes_kept)

    if(n_genes_kept == 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    do i = 1, n_genes
        call logical_as_c_int(mask_logical(i), mask(i))
    end do
end subroutine filter_unassigned_genes_C