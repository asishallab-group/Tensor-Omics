#include <src/macros.h>

!> Parsers for the plain-text input formats TensorOmics data sets are built from (gene-expression
!> TSV/CSV files, OrthoFinder-style family files), plus small array-filtering helpers.
!|
!| These populate the raw `gene_ids` / `expression` / `gene_to_family` / `family_ids` arrays that
!| [[tox_data_archive(module):save_tox_data(subroutine)]] later persists; unlike the archive/serde
!| layers, everything here works from delimited text rather than the library's binary array format.
!|
!| The readers are deliberately simple, not a general CSV reader: a line is read whole up to
!| 4096 characters and refused with `ERR_STRING_TOO_LONG` past that; each delimiter separates two
!| fields, so two in a row enclose an empty one; blanks around a field are not part of it; there
!| is no quoting. Whatever a reader cannot read as it is meant, it refuses rather than guesses.
!| After an error, the contents of every output are undefined.
module tox_data_tools
    use f42_safeguard
    use iso_fortran_env, only: real64, int32, int64, iostat_end, iostat_eor
    use iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, is_err, is_ok, check_io_stat
    use tox_errors, only: ERR_INVALID_INPUT, ERR_FILE_OPEN, ERR_READ_DATA, ERR_SIZE_MISMATCH, ERR_STRING_TOO_LONG
    use tox_errors, only: ERR_EMPTY_INPUT
    use f42_config, only: DEBUG
    M_IMPLICIT_NONE
    private

    public :: read_gene_ids_from_tsv_file
    public :: read_expression_vectors_tsv
    public :: read_orthofinder_file
    public :: split_string
    public :: get_unassigned_mask, apply_unassigned_mask

    integer(int32), parameter :: MAX_LINE_LENGTH = 4096
        !! The longest line the readers take; a longer one is refused with ERR_STRING_TOO_LONG.
    character(len=1), parameter :: TAB = char(9)
    character(len=1), parameter :: CARRIAGE_RETURN = char(13)

contains

    ! ---- reading lines and fields, shared by the readers ----

    !> AUTHOR_AARON_SCHROEDER
    !| Opens `filename` to read it, setting ERR_FILE_OPEN for argument `arg_pos` if it cannot. The
    !| unit is -1 when nothing was opened, so that nothing is closed then.
    subroutine open_text_file(filename, arg_pos, unit, ierr)
        character(len=*), intent(in) :: filename
            !! The file
        integer(int32), intent(in) :: arg_pos
            !! The position of the argument naming the file, for the error
        integer(int32), intent(out) :: unit
            !! The unit, or -1
        integer(int32), intent(inout) :: ierr
            !! Error code
        integer(int32) :: ios

        unit = -1
        open (newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            unit = -1
            call set_err_once(ierr, ERR_FILE_OPEN, arg_pos)
            if (DEBUG) write (*, *) 'Error opening file: ', trim(filename)
        end if
    end subroutine open_text_file

    !> AUTHOR_AARON_SCHROEDER
    !| Skips `n_lines` lines, of any length; a file with fewer is ERR_READ_DATA. A last line without
    !| a line ending counts (an advancing read reports the end of the file for it on gfortran).
    subroutine skip_lines(unit, n_lines, arg_pos, ierr)
        integer(int32), intent(in) :: unit, n_lines, arg_pos
            !! The unit, how many lines, the position of the file argument
        integer(int32), intent(inout) :: ierr
            !! Error code
        integer(int32) :: i_line, ios, n_read
        character(len=256) :: chunk
        logical :: any_read

        do i_line = 1, n_lines
            any_read = .false.
            ! in chunks, without advancing, until the end of the line: so a line of any length is
            ! skipped, and one that has characters counts even where the file ends right after it
            do
                n_read = 0
                read (unit, '(A)', advance='no', size=n_read, iostat=ios) chunk
                if (n_read > 0) any_read = .true.
                if (ios /= 0) exit
            end do
            if (ios == iostat_eor .or. (ios == iostat_end .and. any_read)) cycle
            call set_err_once(ierr, ERR_READ_DATA, arg_pos)
            return
        end do
    end subroutine skip_lines

    !> AUTHOR_AARON_SCHROEDER
    !| Reads the next line into `line(1:line_length)`, without a carriage return that ends it.
    !| `at_end` at the end of the file; a line longer than MAX_LINE_LENGTH is ERR_STRING_TOO_LONG.
    subroutine read_line(unit, line, line_length, at_end, arg_pos, ierr)
        integer(int32), intent(in) :: unit
            !! The unit
        character(len=MAX_LINE_LENGTH + 1), intent(out) :: line
            !! The line, one character longer than any line taken, so a longer one fills it
        integer(int32), intent(out) :: line_length
            !! The length of the line
        logical, intent(out) :: at_end
            !! Whether the file ended before a line
        integer(int32), intent(in) :: arg_pos
            !! The position of the file argument, for an error
        integer(int32), intent(inout) :: ierr
            !! Error code
        integer(int32) :: ios

        at_end = .false.
        line_length = 0
        read (unit, '(A)', advance='no', size=line_length, iostat=ios) line
        if (ios == iostat_end .and. line_length == 0) then
            at_end = .true.
            return
        end if
        if (ios == 0) then
            ! the whole buffer filled without reaching the end of the line
            call set_err_once(ierr, ERR_STRING_TOO_LONG, arg_pos)
            if (DEBUG) write (*, *) 'Error: a line longer than ', MAX_LINE_LENGTH, ' characters'
            return
        end if
        if (ios /= iostat_eor .and. ios /= iostat_end) then
            call set_err_once(ierr, ERR_READ_DATA, arg_pos)
            return
        end if
        ! dead on gfortran, which drops the carriage return of a CR LF line and of a lone CR ending
        ! the file itself, and on ifx for CR LF; kept for a compiler that does not
        if (line_length > 0) then
            if (line(line_length:line_length) == CARRIAGE_RETURN) line_length = line_length - 1
        end if
    end subroutine read_line

    !> AUTHOR_AARON_SCHROEDER
    !| Reads the next line that is not blank. With `skip_blank`, blank lines are passed over
    !| wherever they are; without, only at the end of the file, and a blank line followed by one
    !| that is not is ERR_INVALID_INPUT (a row of no field). `line_number` counts the lines read.
    subroutine read_data_line(unit, skip_blank, line, line_length, line_number, at_end, arg_pos, ierr)
        integer(int32), intent(in) :: unit
            !! The unit
        logical, intent(in) :: skip_blank
            !! Whether a blank line between two others is passed over
        character(len=MAX_LINE_LENGTH + 1), intent(out) :: line
            !! The line
        integer(int32), intent(out) :: line_length
            !! The length of the line
        integer(int32), intent(inout) :: line_number
            !! The number of the last line read
        logical, intent(out) :: at_end
            !! Whether the file ended before a line that is not blank
        integer(int32), intent(in) :: arg_pos
            !! The position of the file argument, for an error
        integer(int32), intent(inout) :: ierr
            !! Error code
        integer(int32) :: first_blank

        first_blank = 0
        do
            call read_line(unit, line, line_length, at_end, arg_pos, ierr)
            if (is_err(ierr) .or. at_end) return
            line_number = line_number + 1
            ! a line of blanks and tabs is blank too, as spreadsheets write them
            if (verify(line(1:line_length), ' '//TAB) > 0) exit
            if (first_blank == 0) first_blank = line_number
        end do
        if (first_blank > 0 .and. .not. skip_blank) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
            if (DEBUG) write (*, *) 'Error: line ', first_blank, ' is blank'
        end if
    end subroutine read_data_line

    !> AUTHOR_AARON_SCHROEDER
    !| The fields of `line` between `delimiter`s, as their first and last positions, blanks around
    !| them excluded unless `keep_blanks`: each delimiter separates two fields, so a field may be
    !| empty (last < first), and a line has one field more than delimiters.
    pure subroutine split_fields(line, delimiter, firsts, lasts, n_fields, keep_blanks)
        character(len=*), intent(in) :: line
            !! The line
        character(len=1), intent(in) :: delimiter
            !! The delimiter
        integer(int32), allocatable, intent(out) :: firsts(:), lasts(:)
            !! The first and last position of each field
        integer(int32), intent(out) :: n_fields
            !! How many fields
        logical, intent(in), optional :: keep_blanks
            !! Whether blanks around a field are part of it; they are not by default
        integer(int32) :: i, first
        logical :: trim_blanks

        trim_blanks = .true.
        if (present(keep_blanks)) trim_blanks = .not. keep_blanks
        n_fields = 1
        do i = 1, len(line)
            if (line(i:i) == delimiter) n_fields = n_fields + 1
        end do
        allocate (firsts(n_fields), lasts(n_fields))
        n_fields = 0
        first = 1
        do i = 1, len(line) + 1
            if (i <= len(line)) then
                if (line(i:i) /= delimiter) cycle
            end if
            n_fields = n_fields + 1
            firsts(n_fields) = first
            lasts(n_fields) = i - 1
            first = i + 1
        end do
        if (.not. trim_blanks) return
        do i = 1, n_fields
            do while (firsts(i) <= lasts(i))
                if (line(firsts(i):firsts(i)) /= ' ') exit
                firsts(i) = firsts(i) + 1
            end do
            do while (lasts(i) >= firsts(i))
                if (line(lasts(i):lasts(i)) /= ' ') exit
                lasts(i) = lasts(i) - 1
            end do
        end do
    end subroutine split_fields

    !> AUTHOR_AARON_SCHROEDER
    !| Reads `text` as one complete real number: refused (`ios` not 0) when it holds a character
    !| at which a list-directed read would stop early or that it reads specially -- a blank, tab,
    !| comma, semicolon, slash or asterisk --, a sign other than at the start or right after an
    !| exponent letter (a list-directed read takes 2-3 for 2e-3), an exponent letter not followed
    !| by a digit or a sign and a digit (ifx reads 1e as 1), or when the read fails.
    subroutine read_number(text, value, ios)
        character(len=*), intent(in) :: text
            !! The field
        real(real64), intent(out) :: value
            !! Its value
        integer(int32), intent(out) :: ios
            !! 0 if it was read
        character(len=*), parameter :: EXPONENT_LETTERS = 'eEdD', SIGNS = '+-', DIGITS = '0123456789'
        integer(int32) :: i, digit_at

        value = 0
        ios = 1
        if (len(text) == 0 .or. scan(text, ' ,;/*'//TAB) > 0) return
        do i = 2, len(text)
            if (index(SIGNS, text(i:i)) > 0 .and. index(EXPONENT_LETTERS, text(i - 1:i - 1)) == 0) return
        end do
        do i = 1, len(text)
            if (index(EXPONENT_LETTERS, text(i:i)) == 0) cycle
            digit_at = i + 1
            if (digit_at <= len(text)) then
                if (index(SIGNS, text(digit_at:digit_at)) > 0) digit_at = digit_at + 1
            end if
            if (digit_at > len(text)) return
            if (index(DIGITS, text(digit_at:digit_at)) == 0) return
        end do
        read (text, *, iostat=ios) value
    end subroutine read_number

    ! ---- the readers ----

    !> M_EXPORT_C
    !| summary: Read expression vectors from csv/tsv files
    !| AUTHOR_AARON_SCHROEDER
    !|
    !| Each file's rows give the values of the gene their `gene_col` names: the `value_cols` of a
    !| row go to rows `start_row`, `start_row + 1`, ... of that gene's column, and each further file
    !| continues below the last. A row naming a gene not in `gene_ids` is skipped; a gene may be
    !| named by one row of a file. Every row must have every column read, none of them empty; a
    !| value must be one complete, finite number. Blank lines at the end of a file are ignored; a
    !| file of no data row is refused with ERR_EMPTY_INPUT. A line may be 4096 characters long.
    !| Values of genes no row names are left as they were; after an error, `expression_vectors`
    !| is undefined.
    subroutine read_expression_vectors_tsv(file_list, gene_ids, expression_vectors, &
                                           n_header_rows, gene_col, value_cols, start_row, delimiter, ierr)
        use f42_xxh3_hashmap
        use ieee_arithmetic, only: ieee_is_finite
        character(len=*), intent(in) :: file_list(:)
            !! List of files to read from
        character(len=*), intent(in) :: gene_ids(:)
            !! Array of gene IDS
        real(real64), intent(inout) :: expression_vectors(:, :)
            !! Array of expression vectors
        integer(int32), intent(in) :: n_header_rows
            !! Number of header rows to skip
        integer(int32), intent(in) :: gene_col
            !! Index of column with gene_ids
        integer(int32), intent(in) :: value_cols(:)
            !! Indicies of columns containing values
        integer(int32), intent(in) :: start_row
            !! Row in the expression vectors to start in
        character(len=1), intent(in), optional :: delimiter
            !! optional delimiter
            !! DM_DEFAULT(char(9))
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_file, i_gene, k, unit, ios, idx, n_genes, n_value_cols, n_fields, n_columns
        integer(int32) :: current_sample, line_length, line_number, n_rows
        integer(int32), allocatable :: firsts(:), lasts(:)
        logical, allocatable :: named(:)
        character(len=MAX_LINE_LENGTH + 1) :: line
        character(len=1) :: actual_delimiter
        real(real64) :: value
        logical :: at_end
        type(hashmap_type) :: gene_map

        call set_ok(ierr)
        n_genes = size(gene_ids)
        n_value_cols = size(value_cols)
        M_DEFAULT_VAL(delimiter, actual_delimiter, TAB)

        ! the argument positions name the offending argument, as the generated wrappers do
        if (n_value_cols <= 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT, 6)
            if (DEBUG) write (*, *) 'Error: No value columns specified.'
            return
        end if
        if (any(value_cols < 1)) then
            call set_err_once(ierr, ERR_INVALID_INPUT, 6)
            return
        end if
        if (gene_col < 1) then
            call set_err_once(ierr, ERR_INVALID_INPUT, 5)
            return
        end if
        if (start_row < 1) then
            call set_err_once(ierr, ERR_INVALID_INPUT, 7)
            return
        end if
        if (size(expression_vectors, 2) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 3)
            if (DEBUG) write (*, *) 'Error: expression_vectors column count does not match n_genes.'
            return
        end if
        ! every file's value columns go below the previous file's, from start_row on
        if (int(start_row, int64) - 1 + int(size(file_list), int64)*n_value_cols > size(expression_vectors, 1)) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 3)
            return
        end if
        n_columns = max(gene_col, maxval(value_cols))

        allocate (named(n_genes), stat=ios)
        call check_io_stat(ios, ierr)
        if (is_err(ierr)) return
        call hashmap_create(gene_map, initial_size=n_genes)
        do i_gene = 1, n_genes
            call hashmap_put(gene_map, trim(gene_ids(i_gene)), i_gene)
        end do

        current_sample = start_row - 1
        unit = -1
        files: do i_file = 1, size(file_list)
            if (DEBUG) write (*, *) 'Reading file: ', trim(file_list(i_file))
            call open_text_file(file_list(i_file), 1, unit, ierr)
            if (is_err(ierr)) exit files
            call skip_lines(unit, n_header_rows, 1, ierr)
            line_number = max(n_header_rows, 0)
            n_rows = 0
            named = .false.
            rows: do while (is_ok(ierr))
                call read_data_line(unit, .false., line, line_length, line_number, at_end, 1, ierr)
                if (is_err(ierr) .or. at_end) exit rows
                n_rows = n_rows + 1
                call split_fields(line(1:line_length), actual_delimiter, firsts, lasts, n_fields)
                if (n_fields < n_columns) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                    if (DEBUG) write (*, *) 'Error: line ', line_number, ' has ', n_fields, ' fields'
                    exit rows
                end if
                if (lasts(gene_col) < firsts(gene_col)) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                    if (DEBUG) write (*, *) 'Error: line ', line_number, ' has an empty gene id'
                    exit rows
                end if
                idx = hashmap_get(gene_map, line(firsts(gene_col):lasts(gene_col)))
                if (idx == -1) then
                    if (DEBUG) write (*, *) 'Warning: Gene ', line(firsts(gene_col):lasts(gene_col)), &
                        ' not found in master gene list'
                    cycle rows
                end if
                if (named(idx)) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                    if (DEBUG) write (*, *) 'Error: line ', line_number, ' names its gene a second time'
                    exit rows
                end if
                named(idx) = .true.
                do k = 1, n_value_cols
                    associate (first => firsts(value_cols(k)), last => lasts(value_cols(k)))
                        if (last < first) then
                            call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                            if (DEBUG) write (*, *) 'Error: line ', line_number, ' column ', value_cols(k), ' is empty'
                            exit rows
                        end if
                        call read_number(line(first:last), value, ios)
                    end associate
                    if (ios /= 0) then
                        call set_err_once(ierr, ERR_READ_DATA, 1)
                        if (DEBUG) write (*, *) 'Error: line ', line_number, ' column ', value_cols(k), ' is no number'
                        exit rows
                    end if
                    if (.not. ieee_is_finite(value)) then
                        call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                        if (DEBUG) write (*, *) 'Error: line ', line_number, ' column ', value_cols(k), ' is not finite'
                        exit rows
                    end if
                    expression_vectors(current_sample + k, idx) = value
                end do
            end do rows
            if (is_ok(ierr) .and. n_rows == 0) call set_err_once(ierr, ERR_EMPTY_INPUT, 1)
            close (unit)
            unit = -1
            if (is_err(ierr)) exit files
            current_sample = current_sample + n_value_cols
        end do files

        call hashmap_destroy(gene_map)
    end subroutine read_expression_vectors_tsv

    !> M_EXPORT_C
    !| summary: Only read the gene ids from a tsv file
    !| AUTHOR_AARON_SCHROEDER
    !|
    !| Each tab-separated row after the header gives one id, from its `gene_col`, in file order.
    !| A row without that column, or with it empty, is refused, and so is an id longer than the
    !| strings of `gene_ids`. Blank lines at the end of the file are ignored; a file of no data row
    !| is refused with ERR_EMPTY_INPUT. A line may be 4096 characters long. With fewer rows than
    !| `gene_ids` has elements, the ones past the last row are left undefined, and the number of
    !| rows read is not returned: size `gene_ids` to the file. After an error, `gene_ids` is
    !| undefined.
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

        integer(int32) :: unit, row_count, line_length, line_number, n_fields
        integer(int32), allocatable :: firsts(:), lasts(:)
        character(len=MAX_LINE_LENGTH + 1) :: line
        logical :: at_end

        call set_ok(ierr)
        row_count = 0

        if (gene_col < 1) then
            call set_err_once(ierr, ERR_INVALID_INPUT, 4)
            return
        end if

        call open_text_file(filename, 1, unit, ierr)
        if (is_err(ierr)) return
        call skip_lines(unit, n_header_rows, 1, ierr)
        line_number = max(n_header_rows, 0)
        do while (is_ok(ierr))
            call read_data_line(unit, .false., line, line_length, line_number, at_end, 1, ierr)
            if (is_err(ierr) .or. at_end) exit
            call split_fields(line(1:line_length), TAB, firsts, lasts, n_fields)
            if (n_fields < gene_col) then
                call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                if (DEBUG) write (*, *) 'Error: line ', line_number, ' has no column ', gene_col
                exit
            end if
            if (lasts(gene_col) < firsts(gene_col)) then
                call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                if (DEBUG) write (*, *) 'Error: line ', line_number, ' has an empty gene id'
                exit
            end if
            if (lasts(gene_col) - firsts(gene_col) + 1 > len(gene_ids)) then
                call set_err_once(ierr, ERR_STRING_TOO_LONG, 2)
                if (DEBUG) write (*, *) 'Error: line ', line_number, ' has an id longer than ', len(gene_ids)
                exit
            end if
            row_count = row_count + 1
            if (row_count > size(gene_ids)) then
                call set_err_once(ierr, ERR_SIZE_MISMATCH, 2)  ! More genes than allocated space
                exit
            end if
            gene_ids(row_count) = line(firsts(gene_col):lasts(gene_col))
        end do
        close (unit)
        if (is_ok(ierr) .and. row_count == 0) call set_err_once(ierr, ERR_EMPTY_INPUT, 1)
    end subroutine read_gene_ids_from_tsv_file

    !> M_EXPORT_C
    !| summary: Read a family file (Orthofinder)
    !| AUTHOR_AARON_SCHROEDER
    !|
    !| After its header line, each line of an OrthoFinder `Orthogroups.tsv` is a family: its id,
    !| then a tab-separated column per species, each a comma-separated list of genes (empty
    !| columns allowed). Family k is the k-th such line, and each gene of `gene_ids` gets the
    !| index of the family that lists it, or M_GENE_TO_FAM_SENTINEL when none does; listed genes
    !| not in `gene_ids` are ignored. Blank lines are skipped. A family id longer than the strings
    !| of `family_ids`, a line without one, and more families than `family_ids` holds are refused.
    !| A line may be 4096 characters long. With fewer families than `family_ids` has elements, the
    !| ones past the last family are left undefined; after an error, all outputs are undefined.
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

        integer(int32) :: unit, i_gene, i_field, i_listed, fam_idx, gene_idx, line_length, line_number
        integer(int32) :: n_fields, n_listed
        integer(int32), allocatable :: firsts(:), lasts(:), listed_firsts(:), listed_lasts(:)
        character(len=MAX_LINE_LENGTH + 1) :: line
        logical :: at_end
        type(hashmap_type) :: gene_map

        call set_ok(ierr)
        if (size(gene_to_fam) /= size(gene_ids)) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 4)
            return
        end if
        gene_to_fam = 0

        call open_text_file(filename, 1, unit, ierr)
        if (is_err(ierr)) return

        call hashmap_create(gene_map, initial_size=size(gene_ids))
        do i_gene = 1, size(gene_ids)
            call hashmap_put(gene_map, gene_ids(i_gene), i_gene)
        end do

        call skip_lines(unit, 1, 1, ierr)
        line_number = 1
        fam_idx = 0
        do while (is_ok(ierr))
            call read_data_line(unit, .true., line, line_length, line_number, at_end, 1, ierr)
            if (is_err(ierr) .or. at_end) exit
            call split_fields(line(1:line_length), TAB, firsts, lasts, n_fields)
            if (lasts(1) < firsts(1)) then
                call set_err_once(ierr, ERR_INVALID_INPUT, 1)
                if (DEBUG) write (*, *) 'Error: line ', line_number, ' has no family id'
                exit
            end if
            if (lasts(1) - firsts(1) + 1 > len(family_ids)) then
                call set_err_once(ierr, ERR_STRING_TOO_LONG, 3)
                exit
            end if
            fam_idx = fam_idx + 1
            if (fam_idx > size(family_ids)) then
                call set_err_once(ierr, ERR_SIZE_MISMATCH, 3)
                exit
            end if
            family_ids(fam_idx) = line(firsts(1):lasts(1))

            do i_field = 2, n_fields
                if (lasts(i_field) < firsts(i_field)) cycle
                associate (column => line(firsts(i_field):lasts(i_field)))
                    call split_fields(column, ',', listed_firsts, listed_lasts, n_listed)
                    do i_listed = 1, n_listed
                        if (listed_lasts(i_listed) < listed_firsts(i_listed)) cycle
                        gene_idx = hashmap_get(gene_map, column(listed_firsts(i_listed):listed_lasts(i_listed)))
                        if (gene_idx > 0) then
                            gene_to_fam(gene_idx) = fam_idx
                        else
                            if (DEBUG) write (*, *) 'Gene not found in hashmap: ', &
                                column(listed_firsts(i_listed):listed_lasts(i_listed))
                        end if
                    end do
                end associate
            end do
        end do
        close (unit)
        call hashmap_destroy(gene_map)
    end subroutine read_orthofinder_file

    !> AUTHOR_AARON_SCHROEDER
    !| Helper subroutine to split strings: the fields of `input` between delimiters, each without
    !| the blanks around it; a run of delimiters separates once, and delimiters at the start or
    !| the end separate nothing, so no field is empty unless it held only blanks.
    subroutine split_string(input, output, delimiter, ierr)
        character(len=*), intent(in) :: input
        !! Input string to split
        character(len=:), allocatable, intent(out) :: output(:)
        !! Output of parts
        character(len=1), optional, intent(in) :: delimiter
        !! delimiter, default is space
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32) :: ios, n, i, n_fields, raw_first, raw_last
        integer(int32), allocatable :: firsts(:), lasts(:), kept(:)
        character(len=1) :: delim

        call set_ok(ierr)
        M_DEFAULT_VAL(delimiter, delim, ' ')

        ! the raw fields, blanks included, of the input without its trailing blanks; one is kept
        ! when it has any character at all, as a run of delimiters separates once
        call split_fields(input(1:len_trim(input)), delim, firsts, lasts, n_fields, keep_blanks=.true.)
        allocate (kept(n_fields), stat=ios)
        call check_io_stat(ios, ierr)
        if (is_err(ierr)) return
        n = 0
        do i = 1, n_fields
            if (lasts(i) >= firsts(i)) then
                n = n + 1
                kept(n) = i
            end if
        end do

        if (n == 0) then
            allocate (character(len=0) :: output(0), stat=ios)
        else
            allocate (character(len=maxval(lasts(kept(1:n)) - firsts(kept(1:n)) + 1)) :: output(n), stat=ios)
        end if
        call check_io_stat(ios, ierr)
        if (is_err(ierr)) return
        do i = 1, n
            raw_first = firsts(kept(i))
            raw_last = lasts(kept(i))
            output(i) = trim(adjustl(input(raw_first:raw_last)))
        end do
    end subroutine split_string

    !> M_EXPORT_C
    !| summary: Helper to create a mask of genes that are unassigned
    !| AUTHOR_AARON_SCHROEDER
    subroutine get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
        use iso_fortran_env, only: int32
        M_IMPLICIT_NONE

        integer(int32), intent(in) :: gene_to_fam(:)
        !! gene to family mapping
        logical(c_bool), intent(out) :: mask(size(gene_to_fam))
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

    !> AUTHOR_AARON_SCHROEDER
    !| Applies a precomputed mask to filter out unassigned genes
    subroutine apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)

        character(len=*), allocatable, intent(inout) :: gene_ids(:)
            !! gene ids array
        real(real64), allocatable, intent(inout) :: expression_vectors(:, :)
            !! expression vectors array
        integer(int32), allocatable, intent(inout) :: gene_to_fam(:)
            !! gene to family mapping
        logical(c_bool), intent(in) :: mask(:)
            !! mask for mapping
        integer(int32), intent(out) :: n_genes_kept
            !! number of genes kept
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i, j, n_genes_total, n_samples, ios, ios2, ios3
        character(len=len(gene_ids)), allocatable :: temp_gene_ids(:)
        real(real64), allocatable :: temp_expression_vectors(:, :)
        integer(int32), allocatable :: temp_gene_to_fam(:)

        call set_ok(ierr)
        n_genes_total = size(gene_to_fam)
        n_samples = size(expression_vectors, 1)
        n_genes_kept = 0
        if (size(gene_ids) /= n_genes_total) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 1)
            return
        end if
        if (size(expression_vectors, 2) /= n_genes_total) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 2)
            return
        end if
        if (size(mask) /= n_genes_total) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH, 4)
            return
        end if
        n_genes_kept = count(mask)

        if (n_genes_kept == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        allocate (temp_gene_ids(n_genes_kept), stat=ios)
        call check_io_stat(ios, ierr)
        allocate (temp_expression_vectors(n_samples, n_genes_kept), stat=ios2)
        call check_io_stat(ios2, ierr)
        allocate (temp_gene_to_fam(n_genes_kept), stat=ios3)
        call check_io_stat(ios3, ierr)
        if (is_err(ierr)) return

        j = 1
        do i = 1, n_genes_total
            if (mask(i)) then
                temp_gene_ids(j) = gene_ids(i)
                temp_expression_vectors(:, j) = expression_vectors(:, i)
                temp_gene_to_fam(j) = gene_to_fam(i)
                j = j + 1
            end if
        end do

        deallocate (gene_ids, expression_vectors, gene_to_fam)
        call move_alloc(temp_gene_ids, gene_ids)
        call move_alloc(temp_expression_vectors, expression_vectors)
        call move_alloc(temp_gene_to_fam, gene_to_fam)
    end subroutine

end module tox_data_tools
