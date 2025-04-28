program main
    use BandedTable_mod
    use iso_fortran_env, only: error_unit
    implicit none

    type(Table) :: my_table
    character(len=256) :: filename
    integer :: ierror
    type(TableColumn), pointer :: col_ptr

    ! --- Example 1: Tab-separated file with header ---
    print *, "--- Example 1: data_header.tsv ---"
    ! Create a dummy data file (data_header.tsv)
    open(10, file='data_header.tsv', status='replace', action='write')
    write(10,'(a)') 'ID'//achar(9)//'Name'//achar(9)//'Value'//achar(9)//'Category'
    write(10,'(a)') '1'//achar(9)//'Apple'//achar(9)//'10.5'//achar(9)//'Fruit'
    write(10,'(a)') '2'//achar(9)//'Banana'//achar(9)//'5.2'//achar(9)//'Fruit'
    write(10,'(a)') '3'//achar(9)//'Carrot'//achar(9)//'2.0'//achar(9)//'Vegetable'
    write(10,'(a)') '4'//achar(9)//'Date'//achar(9)//'2.0'//achar(9)//'Fruit' ! Empty value
    write(10,'(a)') '5'//achar(9)//'Eggplant'//achar(9)//'7.0'//achar(9)//'Vegetable' ! Int value in real col
    close(10)

    filename = 'data_header.tsv'

    print*,'generated Table'
    !   call parse_table(filename, my_table, header=.true., separator=achar(9), ierror=ierror)
    print*,'parsed_Table'
    if (ierror == 0) then
        print *, "Table parsed successfully:"
        call print_preview(my_table)

        !Example Lookup
        col_ptr => get_column_by_name(my_table, 'Value')
        if (associated(col_ptr)) then
           print *, "Column 'Value' found: index =", col_ptr%column_index, ", type =", col_ptr%type_code, &
                   ", start_index =", col_ptr%start_index
        else
           print *, "Column 'Value' not found."
        end if

        call destroy_table(my_table) ! Clean up memory
    else
        write(error_unit, '(a, i0)') "Failed to parse table, error code: ", ierror
    end if
    print *, ""


    ! --- Example 2: Comma-separated file without header ---
    print *, "--- Example 2: data_noheader.csv ---"
    ! Create another dummy data file (data_noheader.csv)
    open(10, file='data_noheader.csv', status='replace', action='write')
    write(10,'(a)') '101,SensorA,99.1'
    write(10,'(a)') '102,SensorB,103.5'
    write(10,'(a)') '103,SensorC,-5.0'
    close(10)

    filename = 'data_noheader.csv'
    call parse_table(filename, my_table, header=.false., separator=',', ierror=ierror)

    if (ierror == 0) then
        print *, "Table parsed successfully:"
        call print_preview(my_table)
        call destroy_table(my_table)
    else
        write(error_unit, '(a, i0)') "Failed to parse table, error code: ", ierror
    end if
    print *, ""

end program main