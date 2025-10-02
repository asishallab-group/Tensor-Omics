! filepath: test/mod_test_csv_file_reader.f90
!> Unit test suite for CSV file reader module.
module mod_test_csv_file_reader
  use asserts
  use tox_csv_file_reader
  use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_FILE_EMPTY, ERR_DIM_MISMATCH, is_err, is_ok
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

! // TODO: Check again with real big csv files
! // TODO: Check with "\t" for tab separator

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(14)
    all_tests(1) = test_case("test_all_data_types_with_header", test_all_data_types_with_header)
    all_tests(2) = test_case("test_all_data_types_without_header", test_all_data_types_without_header)
    all_tests(3) = test_case("test_all_data_types_with_column_names", test_all_data_types_with_column_names)
    all_tests(4) = test_case("test_column_names_override_header", test_column_names_override_header)
    all_tests(5) = test_case("test_csv_with_comments", test_csv_with_comments)
    all_tests(6) = test_case("test_tsv_with_tab_separator", test_tsv_with_tab_separator)
    all_tests(7) = test_case("test_inconsistent_separators_error", test_inconsistent_separators_error)
    all_tests(8) = test_case("test_empty_fields_error", test_empty_fields_error)
    all_tests(9) = test_case("test_empty_column_names", test_empty_column_names)
    all_tests(10) = test_case("test_empty_column_types", test_empty_column_types)
    all_tests(11) = test_case("test_empty_csv_file", test_empty_csv_file)
    all_tests(12) = test_case("test_different_line_breaks", test_different_line_breaks)
    all_tests(13) = test_case("test_getter_functions", test_getter_functions)
    all_tests(14) = test_case("test_serialization_deserialization", test_serialization_deserialization)
  end function get_all_tests

  !> Run all CSV file reader tests.
  subroutine run_all_tests_csv_file_reader()
    type(test_case) :: all_tests(14)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All CSV file reader tests passed successfully."
  end subroutine run_all_tests_csv_file_reader

  !> Run specific CSV file reader tests by name.
  subroutine run_named_tests_csv_file_reader(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(14)
    integer(int32) :: i, j
    logical :: found
    all_tests = get_all_tests()
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_csv_file_reader

  !> Test CSV reader with all 5 data types on a 20x10 mixed CSV file with header
  subroutine test_all_data_types_with_header()
    ! Define arrays for all 5 data types - 2 columns of each type (20 rows x 2 cols)
    integer(int32) :: int_cols(20, 2), expected_int_cols(20, 2)
    real(real64) :: real_cols(20, 2), expected_real_cols(20, 2)  
    character(len=64) :: char_cols(20, 2), expected_char_cols(20, 2)
    logical :: logical_cols(20, 2), expected_logical_cols(20, 2)
    complex(real64) :: complex_cols(20, 2), expected_complex_cols(20, 2)
    character(len=64) :: header(10), expected_header(10)
    integer(int32) :: metadata(2, 10), expected_metadata(2, 10)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_reader_temp.csv"
    
    ! Column types: 2 int, 2 real, 2 char, 2 logical, 2 complex
    ! Column order: int, real, char, logical, complex, int, real, char, logical, complex
    integer(int32) :: column_types(10) = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
    
    ! Create temporary test CSV file
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV header
    write(file_unit, '(A)') "IntCol1,RealCol1,CharCol1,LogicalCol1,ComplexCol1,IntCol2,RealCol2,CharCol2,LogicalCol2,ComplexCol2"
    
    ! Write CSV data rows
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0),11,11.5,AA,T,(11.0,12.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0),12,12.5,BB,F,(12.0,13.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0),13,13.5,CC,T,(13.0,14.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0),14,14.5,DD,F,(14.0,15.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0),15,15.5,EE,T,(15.0,16.0)"
    write(file_unit, '(A)') "6,6.5,F,F,(6.0,7.0),16,16.5,FF,F,(16.0,17.0)"
    write(file_unit, '(A)') "7,7.5,G,T,(7.0,8.0),17,17.5,GG,T,(17.0,18.0)"
    write(file_unit, '(A)') "8,8.5,H,F,(8.0,9.0),18,18.5,HH,F,(18.0,19.0)"
    write(file_unit, '(A)') "9,9.5,I,T,(9.0,10.0),19,19.5,II,T,(19.0,20.0)"
    write(file_unit, '(A)') "10,10.5,J,F,(10.0,11.0),20,20.5,JJ,F,(20.0,21.0)"
    write(file_unit, '(A)') "11,11.5,K,T,(11.0,12.0),21,21.5,KK,T,(21.0,22.0)"
    write(file_unit, '(A)') "12,12.5,L,F,(12.0,13.0),22,22.5,LL,F,(22.0,23.0)"
    write(file_unit, '(A)') "13,13.5,M,T,(13.0,14.0),23,23.5,MM,T,(23.0,24.0)"
    write(file_unit, '(A)') "14,14.5,N,F,(14.0,15.0),24,24.5,NN,F,(24.0,25.0)"
    write(file_unit, '(A)') "15,15.5,O,T,(15.0,16.0),25,25.5,OO,T,(25.0,26.0)"
    write(file_unit, '(A)') "16,16.5,P,F,(16.0,17.0),26,26.5,PP,F,(26.0,27.0)"
    write(file_unit, '(A)') "17,17.5,Q,T,(17.0,18.0),27,27.5,QQ,T,(27.0,28.0)"
    write(file_unit, '(A)') "18,18.5,R,F,(18.0,19.0),28,28.5,RR,F,(28.0,29.0)"
    write(file_unit, '(A)') "19,19.5,S,T,(19.0,20.0),29,29.5,SS,T,(29.0,30.0)"
    write(file_unit, '(A)') "20,20.5,T,F,(20.0,21.0),30,30.5,TT,F,(30.0,31.0)"
    
    close(file_unit)

    do i = 1, 20
      expected_int_cols(i, 1) = i         ! Column 1: 1,2,3...20
      expected_int_cols(i, 2) = i + 10    ! Column 6: 11,12,13...30
    end do
    
    do i = 1, 20
      expected_real_cols(i, 1) = real(i, real64) + 0.5_real64  ! 1.5, 2.5, 3.5, ... 20.5
      expected_real_cols(i, 2) = real(i + 10, real64) + 0.5_real64  ! 11.5, 12.5, 13.5, ... 30.5
    end do
    
    expected_char_cols(:, :) = reshape([character(len=64) :: "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
                                                     "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
                                                     "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", &
                                                     "KK", "LL", "MM", "NN", "OO", "PP", "QQ", "RR", "SS", "TT"], [20,2])
    
    expected_logical_cols(:, :) = reshape([.true., .false., .true., .false., .true., .false., .true., .false., .true., .false., &
                                   .true., .false., .true., .false., .true., .false., .true., .false., .true., .false., &
                                   .true., .false., .true., .false., .true., .false., .true., .false., .true., .false., &
                                   .true., .false., .true., .false., .true., .false., .true., .false., .true., .false.], [20,2])
    
    do i = 1, 20
      expected_complex_cols(i, 1) = cmplx(real(i, real64), real(i + 1, real64), real64)        ! (1,2), (2,3), ...
      expected_complex_cols(i, 2) = cmplx(real(i + 10, real64), real(i + 11, real64), real64)  ! (11,12), (12,13), ...
    end do

    expected_header = [character(len=64) :: "IntCol1", "RealCol1", "CharCol1", "LogicalCol1", "ComplexCol1", &
                       "IntCol2", "RealCol2", "CharCol2", "LogicalCol2", "ComplexCol2"]

    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 1, 2, 2, 2, 3, 2, 4, 2, 5, 2], [2, 10])

    ! Read the CSV file
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading was successful and compare all data to expected values
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed")

    call assert_equal_array_int(int_cols, expected_int_cols, 40, "Integer columns should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 40, 1.0e-10_real64, "Real columns should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 40, "Character columns should match")
    call assert_equal_array_int(merge(1, 0, logical_cols), merge(1, 0, expected_logical_cols), 40, "Logical columns should match")
    call assert_equal_array_complex(complex_cols, expected_complex_cols, 40, 1.0e-10_real64, "Complex columns should match")
    call assert_equal_array_char(header, expected_header, 64, 10, "Expected headers should match")
    call assert_equal_array_int(metadata, expected_metadata, 20, "Expected metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_all_data_types_with_header

  !> Test CSV reader with all 5 data types on a 5x5 CSV file without header
  subroutine test_all_data_types_without_header()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(5), expected_header(5)
    integer(int32) :: metadata(2, 5), expected_metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_no_header_temp.csv"
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Create temporary test CSV file (without header)
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV data rows (no header line)
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0)"
    
    close(file_unit)
    
    ! Set up expected values for validation
    do i = 1, 5
      expected_int_cols(i, 1) = i         ! 1,2,3,4,5
    end do
    
    do i = 1, 5
      expected_real_cols(i, 1) = real(i, real64) + 0.5_real64  ! 1.5, 2.5, 3.5, 4.5, 5.5
    end do
    
    expected_char_cols(:, 1) = [character(len=64) :: "A", "B", "C", "D", "E"]
    
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    
    do i = 1, 5
      expected_complex_cols(i, 1) = cmplx(real(i, real64), real(i + 1, real64), real64)  ! (1,2), (2,3), (3,4), (4,5), (5,6)
    end do

    expected_header = [character(len=64) :: "1_int", "2_real", "3_char", "4_logical", "5_complex"]
    
    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1], [2, 5])

    ! Read the CSV file (has_header = .false.)
    call read_table(test_file, column_types, .false., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading was successful
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed")
    
    ! Validate all data columns
    call assert_equal_array_int(int_cols, expected_int_cols, 5, "Integer column should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 5, 1.0e-10_real64, "Real column should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 5, "Character column should match")
    call assert_equal_array_int(merge(1, 0, logical_cols), merge(1, 0, expected_logical_cols), 5, "Logical column should match")
    call assert_equal_array_complex(complex_cols, expected_complex_cols, 5, 1.0e-10_real64, "Complex column should match")
    ! Validate auto-generated headers
    call assert_equal_array_char(header, expected_header, 64, 5, "Auto-generated headers should match")
    
    ! Validate metadata
    call assert_equal_array_int(metadata, expected_metadata, 10, "Expected metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_all_data_types_without_header

  !> Test CSV reader with all 5 data types on a 5x5 CSV file without header but with custom column names
  subroutine test_all_data_types_with_column_names()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(5), expected_header(5)
    integer(int32) :: metadata(2, 5), expected_metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_custom_names_temp.csv"
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Custom column names to override auto-generated ones
    character(len=64), dimension(5) :: custom_column_names = [character(len=64) :: &
        "MyInteger", "MyReal", "MyCharacter", "MyLogical", "MyComplex"]
    
    ! Create temporary test CSV file (without header)
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV data rows (no header line)
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0)"
    
    close(file_unit)
    
    ! Set up expected values for validation
    do i = 1, 5
      expected_int_cols(i, 1) = i         ! 1,2,3,4,5
    end do
    
    do i = 1, 5
      expected_real_cols(i, 1) = real(i, real64) + 0.5_real64  ! 1.5, 2.5, 3.5, 4.5, 5.5
    end do
    
    expected_char_cols(:, 1) = [character(len=64) :: "A", "B", "C", "D", "E"]
    
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    
    do i = 1, 5
      expected_complex_cols(i, 1) = cmplx(real(i, real64), real(i + 1, real64), real64)  ! (1,2), (2,3), (3,4), (4,5), (5,6)
    end do

    expected_header = custom_column_names
    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1], [2, 5])

    ! Read the CSV file (has_header = .false., but with custom column_names)
    call read_table(test_file, column_types, .false., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr, ',', custom_column_names)
    
    ! Check that reading was successful
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed")
    
    ! Validate all data columns
    call assert_equal_array_int(int_cols, expected_int_cols, 5, "Integer column should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 5, 1.0e-10_real64, "Real column should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 5, "Character column should match")
    call assert_equal_array_int(merge(1, 0, logical_cols), merge(1, 0, expected_logical_cols), 5, "Logical column should match")
    call assert_equal_array_complex(complex_cols, expected_complex_cols, 5, 1.0e-10_real64, "Complex column should match")
    
    ! Validate custom headers (should override auto-generated ones)
    call assert_equal_array_char(header, expected_header, 64, 5, "Custom column names should be used")
    ! Validate metadata
    call assert_equal_array_int(metadata, expected_metadata, 10, "Expected metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_all_data_types_with_column_names

  !> Test CSV reader where custom column names override CSV file header
  subroutine test_column_names_override_header()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(5), expected_header(5)
    integer(int32) :: metadata(2, 5), expected_metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_override_header_temp.csv"
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Custom column names that should override the CSV file header
    character(len=64), dimension(5) :: override_column_names = [character(len=64) :: &
        "OverrideInt", "OverrideReal", "OverrideChar", "OverrideLogical", "OverrideComplex"]
    
    ! Create temporary test CSV file (WITH header that should be ignored)
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV header (these names should be ignored in favor of override_column_names)
    write(file_unit, '(A)') "IgnoredInt,IgnoredReal,IgnoredChar,IgnoredLogical,IgnoredComplex"
    
    ! Write CSV data rows
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0)"
    
    close(file_unit)
    
    ! Set up expected values for validation
    do i = 1, 5
      expected_int_cols(i, 1) = i         ! 1,2,3,4,5
    end do
    
    do i = 1, 5
      expected_real_cols(i, 1) = real(i, real64) + 0.5_real64  ! 1.5, 2.5, 3.5, 4.5, 5.5
    end do
    
    expected_char_cols(:, 1) = [character(len=64) :: "A", "B", "C", "D", "E"]
    
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    
    do i = 1, 5
      expected_complex_cols(i, 1) = cmplx(real(i, real64), real(i + 1, real64), real64)  ! (1,2), (2,3), (3,4), (4,5), (5,6)
    end do

    ! Expected headers (should be our override names, NOT the file header)
    expected_header = override_column_names
    ! Validate metadata
    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1], [2, 5])

    ! Read the CSV file (has_header = .true., but column_names should override the file header)
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr, ',', override_column_names)
    
    ! Check that reading was successful
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed")
    
    ! Validate all data columns
    call assert_equal_array_int(int_cols, expected_int_cols, 5, "Integer column should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 5, 1.0e-10_real64, "Real column should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 5, "Character column should match")
    call assert_equal_array_int(merge(1, 0, logical_cols), merge(1, 0, expected_logical_cols), 5, "Logical column should match")
    call assert_equal_array_complex(complex_cols, expected_complex_cols, 5, 1.0e-10_real64, "Complex column should match")
    
    ! Validate that override names are used (NOT the file header names)
    call assert_equal_array_char(header, expected_header, 64, 5, "Override column names should be used instead of file header")
    
    ! Validate metadata
    call assert_equal_array_int(metadata, expected_metadata, 10, "Expected metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_column_names_override_header

  !> Test with comments at the start of the file (lines starting with #)
  subroutine test_csv_with_comments()
    ! Define arrays for mixed data types - 3 columns (5 rows x 3 cols)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 0)  ! No logical columns in this test
    complex(real64) :: complex_cols(5, 0)  ! No complex columns in this test
    character(len=64) :: header(3), expected_header(3)
    integer(int32) :: metadata(2, 3), expected_metadata(2, 3)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_comments_temp.csv"
    
    ! Column types: int, real, char
    integer(int32) :: column_types(3) = [1, 2, 3]
    
    ! Create temporary test CSV file with comments at the start
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file with comments"
      stop 1
    end if
    
    ! Write comment lines at the start
    write(file_unit, '(A)') "# This is a test CSV file with comments"
    write(file_unit, '(A)') "# Author: Test Suite"
    write(file_unit, '(A)') "# Date: 2025-10-02"
    write(file_unit, '(A)') "# Description: Sample data for testing CSV reader with comments"
    write(file_unit, '(A)') "#"
    write(file_unit, '(A)') "# Column descriptions:"
    write(file_unit, '(A)') "# ID: Integer identifier"
    write(file_unit, '(A)') "# Value: Real number value"
    write(file_unit, '(A)') "# Name: Character name"
    write(file_unit, '(A)') "#"
    
    ! Write CSV header
    write(file_unit, '(A)') "ID,Value,Name"
    
    ! Write CSV data rows
    write(file_unit, '(A)') "1,1.1,Alpha"
    write(file_unit, '(A)') "2,2.2,Beta"
    write(file_unit, '(A)') "3,3.3,Gamma"
    write(file_unit, '(A)') "4,4.4,Delta"
    write(file_unit, '(A)') "5,5.5,Epsilon"
    
    close(file_unit)

    ! Set expected values
    do i = 1, 5
      expected_int_cols(i, 1) = i
    end do
    
    do i = 1, 5
      expected_real_cols(i, 1) = real(i, real64) + real(i, real64) / 10.0_real64  ! 1.1, 2.2, 3.3, 4.4, 5.5
    end do
    
    expected_char_cols(:, 1) = [character(len=64) :: "Alpha", "Beta", "Gamma", "Delta", "Epsilon"]

    expected_header = [character(len=64) :: "ID", "Value", "Name"]

    expected_metadata = reshape([1, 1, 2, 1, 3, 1], [2, 3])

    ! Read the CSV file (should skip comment lines automatically)
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading was successful and compare all data to expected values
    call assert_equal_int(ierr, ERR_OK, "CSV reading with comments should succeed")

    call assert_equal_array_int(int_cols, expected_int_cols, 5, "Integer column should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 5, 1.0e-10_real64, "Real column should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 5, "Character column should match")
    call assert_equal_array_char(header, expected_header, 64, 3, "Headers should match")
    call assert_equal_array_int(metadata, expected_metadata, 6, "Metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_csv_with_comments

  !> Test CSV reader with TSV file format using tab separator
  subroutine test_tsv_with_tab_separator()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(5), expected_header(5)
    integer(int32) :: metadata(2, 5), expected_metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_data_temp.tsv"
    character(len=1), parameter :: tab_char = char(9)  ! ASCII code 9 is TAB
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Create temporary TSV file with header
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test TSV file"
      stop 1
    end if
    
    ! Write TSV header (tab-separated)
    write(file_unit, '(A)') "IntColumn" // tab_char // "RealColumn" // tab_char // &
                            "CharColumn" // tab_char // "LogicalColumn" // tab_char // "ComplexColumn"
    
    ! Write TSV data rows (tab-separated)
    write(file_unit, '(A)') "1" // tab_char // "1.5" // tab_char // "A" // tab_char // "T" // tab_char // "(1.0,2.0)"
    write(file_unit, '(A)') "2" // tab_char // "2.5" // tab_char // "B" // tab_char // "F" // tab_char // "(2.0,3.0)"
    write(file_unit, '(A)') "3" // tab_char // "3.5" // tab_char // "C" // tab_char // "T" // tab_char // "(3.0,4.0)"
    write(file_unit, '(A)') "4" // tab_char // "4.5" // tab_char // "D" // tab_char // "F" // tab_char // "(4.0,5.0)"
    write(file_unit, '(A)') "5" // tab_char // "5.5" // tab_char // "E" // tab_char // "T" // tab_char // "(5.0,6.0)"
    
    close(file_unit)
    
    ! Set up expected values for validation
    do i = 1, 5
      expected_int_cols(i, 1) = i         ! 1,2,3,4,5
    end do
    
    do i = 1, 5
      expected_real_cols(i, 1) = real(i, real64) + 0.5_real64  ! 1.5, 2.5, 3.5, 4.5, 5.5
    end do
    
    expected_char_cols(:, 1) = [character(len=64) :: "A", "B", "C", "D", "E"]
    
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    
    do i = 1, 5
      expected_complex_cols(i, 1) = cmplx(real(i, real64), real(i + 1, real64), real64)  ! (1,2), (2,3), (3,4), (4,5), (5,6)
    end do

    expected_header = [character(len=64) :: "IntColumn", "RealColumn", "CharColumn", "LogicalColumn", "ComplexColumn"]
    
    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1], [2, 5])

    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr, tab_char)
    
    ! Check that reading was successful
    call assert_equal_int(ierr, ERR_OK, "TSV reading should succeed")
    
    ! Validate all data columns
    call assert_equal_array_int(int_cols, expected_int_cols, 5, "Integer column should match")
    call assert_equal_array_real(real_cols, expected_real_cols, 5, 1.0e-10_real64, "Real column should match")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 5, "Character column should match")
    call assert_equal_array_int(merge(1, 0, logical_cols), merge(1, 0, expected_logical_cols), 5, "Logical column should match")
    call assert_equal_array_complex(complex_cols, expected_complex_cols, 5, 1.0e-10_real64, "Complex column should match")
    
    ! Validate headers from TSV file
    call assert_equal_array_char(header, expected_header, 64, 5, "TSV header names should be read correctly")
    
    ! Validate metadata
    call assert_equal_array_int(metadata, expected_metadata, 10, "Expected metadata should match")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_tsv_with_tab_separator

  !> Test CSV reader error handling with inconsistent separators
  subroutine test_inconsistent_separators_error()
    integer(int32) :: int_cols(3, 1)
    real(real64) :: real_cols(3, 1)  
    character(len=64) :: char_cols(3, 2)
    logical :: logical_cols(3, 1)
    complex(real64) :: complex_cols(0, 0)
    character(len=64) :: header(4)
    integer(int32) :: metadata(2, 4)
    integer(int32) :: ierr, file_unit, io_status
    character(len=*), parameter :: test_file = "test_inconsistent_separators_temp.csv"
    
    ! Column types: int, real, char, char (4 columns total)
    integer(int32) :: column_types(4) = [1, 2, 3, 3]
    
    ! Create temporary CSV file with inconsistent separators
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV header with consistent separators (comma)
    write(file_unit, '(A)') "IntCol,RealCol,CharCol1,CharCol2"
    
    ! Write data rows with INCONSISTENT separators
    write(file_unit, '(A)') "1|1.5|A,B"
    write(file_unit, '(A)') "2|2.5,C;D"
    write(file_unit, '(A)') "3|3.5|E|F"
    
    close(file_unit)

    ! Try to read the CSV file with comma separator (should fail due to inconsistent separators)
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading failed with an appropriate error code
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "CSV reading should fail due to inconsistent separators")
    
    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_inconsistent_separators_error

  !> Test CSV reader error handling with empty fields
  subroutine test_empty_fields_error()
    ! Define minimal arrays for the test
    integer(int32) :: int_cols(3, 1)
    real(real64) :: real_cols(3, 1)  
    character(len=64) :: char_cols(3, 1)
    logical :: logical_cols(3, 1)
    complex(real64) :: complex_cols(0, 0)  ! No complex columns in this test
    character(len=64) :: header(4)
    integer(int32) :: metadata(2, 4)
    integer(int32) :: ierr, file_unit, io_status
    character(len=*), parameter :: test_file = "test_empty_fields_temp.csv"
    
    ! Column types: int, real, char, logical (4 columns total)
    integer(int32) :: column_types(4) = [1, 2, 3, 4]
    
    ! Create temporary CSV file with empty fields
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV header
    write(file_unit, '(A)') "IntCol,RealCol,CharCol,LogicalCol"
    
    ! Write data rows with EMPTY fields
    write(file_unit, '(A)') ",1.5,A,T"
    write(file_unit, '(A)') "2,,B,F"
    write(file_unit, '(A)') "3,3.5,,T"
    write(file_unit, '(A)') "4,4.5,D,"
    write(file_unit, '(A)') ",,E,"
    
    close(file_unit)

    ! Try to read the CSV file with empty fields (should fail)
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading failed with an appropriate error code
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "CSV reading should fail due to empty fields")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (io_status == 0) then
      close(file_unit, status='delete')
    end if

  end subroutine test_empty_fields_error

  !> Test CSV reader with empty column names array
  subroutine test_empty_column_names()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1)
    real(real64) :: real_cols(5, 1)
    character(len=64) :: char_cols(5, 1)
    logical :: logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1)
    character(len=64) :: header(5)
    integer(int32) :: metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_temp.csv"
    character(len=64) :: column_names(0)
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Create temporary test CSV file
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV data rows
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0)"
    
    close(file_unit)

    ! Read the CSV file (has_header = .true., but column_names should override the file header)
    call read_table(test_file, column_types, .false., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr, ',', column_names)
    
    ! Check that reading failed with an appropriate error code
    call assert_equal_int(ierr, ERR_DIM_MISMATCH, "CSV reading should fail due to empty column names array")
    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_empty_column_names

  !> Test CSV reader with empty column types array
  subroutine test_empty_column_types()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col)
    integer(int32) :: int_cols(5, 1)
    real(real64) :: real_cols(5, 1)
    character(len=64) :: char_cols(5, 1)
    logical :: logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1)
    character(len=64) :: header(5)
    integer(int32) :: metadata(2, 5)
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_csv_temp.csv"
    character(len=64) :: column_names(5)
    ! Empty column types array
    integer(int32) :: column_types(0)
    
    column_names = [character(len=64) :: "int_col", "real_col", "char_col", "logical_col", "complex_col"]
    
    ! Create temporary test CSV file
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV data rows
    write(file_unit, '(A)') "1,1.5,A,T,(1.0,2.0)"
    write(file_unit, '(A)') "2,2.5,B,F,(2.0,3.0)"
    write(file_unit, '(A)') "3,3.5,C,T,(3.0,4.0)"
    write(file_unit, '(A)') "4,4.5,D,F,(4.0,5.0)"
    write(file_unit, '(A)') "5,5.5,E,T,(5.0,6.0)"
    
    close(file_unit)

    ! Read the CSV file (has_header = .true., but column_names should override the file header)
    call read_table(test_file, column_types, .false., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr, ',', column_names)
    
    ! Check that reading failed with an appropriate error code
    call assert_equal_int(ierr, ERR_DIM_MISMATCH, "CSV reading should fail due to empty column types array")
    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_empty_column_types

    !> Test error handling with empty CSV file
  subroutine test_empty_csv_file()
    ! Define minimal arrays for the test
    integer(int32) :: int_cols(5, 1)
    real(real64) :: real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1)
    logical :: logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1)
    character(len=64) :: header(5)
    integer(int32) :: metadata(2, 5)
    integer(int32) :: ierr, file_unit, io_status
    character(len=*), parameter :: test_file = "test_empty_csv_temp.csv"
    
    ! Column types: int, real, char, logical, complex (5 columns expected)
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Create an empty CSV file
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating empty test CSV file"
      stop 1
    end if
    
    ! Close the file immediately without writing anything (creates empty file)
    close(file_unit)

    ! Try to read the empty CSV file (should fail with ERR_FILE_EMPTY)
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading failed with the specific empty file error code
    call assert_equal_int(ierr, ERR_FILE_EMPTY, "CSV reading should fail for empty file")
    
    ! Also test without header
    call read_table(test_file, column_types, .false., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading failed again for empty file
    call assert_equal_int(ierr, ERR_FILE_EMPTY, "CSV reading should fail for empty file (no header)")

    ! Clean up temporary test file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_empty_csv_file

  !> Test CSV reader with different line break types (LF, CR, CR-LF)
  subroutine test_different_line_breaks()
    ! Define arrays for a simple 3x3 CSV test
    integer(int32) :: int_cols(3, 1), expected_int_cols(3, 1)
    real(real64) :: real_cols(3, 1), expected_real_cols(3, 1)  
    character(len=64) :: char_cols(3, 1), expected_char_cols(3, 1)
    logical :: logical_cols(3, 0)  ! No logical columns in this test
    complex(real64) :: complex_cols(3, 0)  ! No complex columns in this test
    character(len=64) :: header(3), expected_header(3)
    integer(int32) :: metadata(2, 3), expected_metadata(2, 3)
    integer(int32) :: ierr, i, file_unit, io_status
    
    ! Column types: int, real, char
    integer(int32) :: column_types(3) = [1, 2, 3]
    
    ! Define line break characters
    character(len=1), parameter :: LF = char(10)      ! Unix/Linux
    character(len=1), parameter :: CR = char(13)      ! Classic Mac
    character(len=2), parameter :: CRLF = char(13) // char(10)  ! Windows
    
    ! Test files for each line break type
    character(len=*), parameter :: test_file_lf = "test_linebreaks_lf_temp.csv"
    character(len=*), parameter :: test_file_cr = "test_linebreaks_cr_temp.csv"
    character(len=*), parameter :: test_file_crlf = "test_linebreaks_crlf_temp.csv"
    
    ! Set up expected data (same for all line break tests)
    expected_int_cols(:, 1) = [10, 20, 30]
    expected_real_cols(:, 1) = [1.1_real64, 2.2_real64, 3.3_real64]
    expected_char_cols(:, 1) = [character(len=64) :: "Alpha", "Beta", "Gamma"]
    expected_header = [character(len=64) :: "ID", "Value", "Name"]
    expected_metadata = reshape([1, 1, 2, 1, 3, 1], [2, 3])

    ! Test 1: LF line breaks (Unix/Linux style)
    open(newunit=file_unit, file=test_file_lf, status='replace', action='write', &
         access='stream', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating LF test CSV file"
      stop 1
    end if
    
    ! Write data with explicit LF line breaks
    write(file_unit) "ID,Value,Name" // LF
    write(file_unit) "10,1.1,Alpha" // LF
    write(file_unit) "20,2.2,Beta" // LF
    write(file_unit) "30,3.3,Gamma" // LF
    close(file_unit)

    ! Read LF file
    call read_table(test_file_lf, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "CSV reading with LF line breaks should succeed")
    call assert_equal_array_int(int_cols, expected_int_cols, 3, "Integer column should match (LF)")
    call assert_equal_array_real(real_cols, expected_real_cols, 3, 1.0e-10_real64, "Real column should match (LF)")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 3, "Character column should match (LF)")
    call assert_equal_array_char(header, expected_header, 64, 3, "Headers should match (LF)")

    ! Test 2: CR line breaks (Classic Mac style)
    open(newunit=file_unit, file=test_file_cr, status='replace', action='write', &
         access='stream', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating CR test CSV file"
      stop 1
    end if
    
    ! Write data with explicit CR line breaks
    write(file_unit) "ID,Value,Name" // CR
    write(file_unit) "10,1.1,Alpha" // CR
    write(file_unit) "20,2.2,Beta" // CR
    write(file_unit) "30,3.3,Gamma" // CR
    close(file_unit)

    ! Read CR file
    call read_table(test_file_cr, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)

    call assert_equal_int(ierr, ERR_OK, "CSV reading with CR line breaks should succeed")
    call assert_equal_array_int(int_cols, expected_int_cols, 3, "Integer column should match (CR)")
    call assert_equal_array_real(real_cols, expected_real_cols, 3, 1.0e-10_real64, "Real column should match (CR)")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 3, "Character column should match (CR)")
    call assert_equal_array_char(header, expected_header, 64, 3, "Headers should match (CR)")

    ! Test 3: CR-LF line breaks (Windows style)
    open(newunit=file_unit, file=test_file_crlf, status='replace', action='write', &
         access='stream', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating CR-LF test CSV file"
      stop 1
    end if
    
    ! Write data with explicit CR-LF line breaks
    write(file_unit) "ID,Value,Name" // CRLF
    write(file_unit) "10,1.1,Alpha" // CRLF
    write(file_unit) "20,2.2,Beta" // CRLF
    write(file_unit) "30,3.3,Gamma" // CRLF
    close(file_unit)

    ! Read CR-LF file
    call read_table(test_file_crlf, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "CSV reading with CR-LF line breaks should succeed")
    call assert_equal_array_int(int_cols, expected_int_cols, 3, "Integer column should match (CR-LF)")
    call assert_equal_array_real(real_cols, expected_real_cols, 3, 1.0e-10_real64, "Real column should match (CR-LF)")
    call assert_equal_array_char(char_cols, expected_char_cols, 64, 3, "Character column should match (CR-LF)")
    call assert_equal_array_char(header, expected_header, 64, 3, "Headers should match (CR-LF)")

    ! Clean up temporary test files
    open(newunit=file_unit, file=test_file_lf, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=test_file_cr, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=test_file_crlf, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_different_line_breaks

  !> Test getter functions for all data types
  subroutine test_getter_functions()
    ! Define arrays for all 5 data types - 1 column of each type (5 rows x 1 col each)
    integer(int32) :: int_cols(5, 1), expected_int_cols(5, 1)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)  
    character(len=64) :: char_cols(5, 1), expected_char_cols(5, 1)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(5), expected_header(5)
    integer(int32) :: metadata(2, 5), expected_metadata(2, 5)
    
    ! Arrays for testing getters
    integer(int32) :: retrieved_int_column(5)
    real(real64) :: retrieved_real_column(5)
    character(len=64) :: retrieved_char_column(5)
    logical :: retrieved_logical_column(5)
    complex(real64) :: retrieved_complex_column(5)
    
    integer(int32) :: ierr, getter_ierr, i, io_status, file_unit
    character(len=*), parameter :: test_file = "test_getter_temp.csv"
    
    ! Column types: int, real, char, logical, complex
    integer(int32) :: column_types(5) = [1, 2, 3, 4, 5]
    
    ! Create temporary test CSV file with header
    file_unit = 11
    open(file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file"
      stop 1
    end if
    
    ! Write CSV header
    write(file_unit, '(A)') "IntColumn,RealColumn,CharColumn,LogicalColumn,ComplexColumn"
    
    ! Write CSV data rows (5 rows)
    write(file_unit, '(A)') "10,1.1,Alpha,T,(1.0,2.0)"
    write(file_unit, '(A)') "20,2.2,Beta,F,(3.0,4.0)"
    write(file_unit, '(A)') "30,3.3,Gamma,T,(5.0,6.0)"
    write(file_unit, '(A)') "40,4.4,Delta,F,(7.0,8.0)"
    write(file_unit, '(A)') "50,5.5,Epsilon,T,(9.0,10.0)"
    
    close(file_unit)

    ! Set up expected data
    expected_int_cols(:, 1) = [10, 20, 30, 40, 50]
    expected_real_cols(:, 1) = [1.1_real64, 2.2_real64, 3.3_real64, 4.4_real64, 5.5_real64]
    expected_char_cols(:, 1) = [character(len=64) :: "Alpha", "Beta", "Gamma", "Delta", "Epsilon"]
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    expected_complex_cols(:, 1) = [cmplx(1.0_real64, 2.0_real64, real64), &
                                   cmplx(3.0_real64, 4.0_real64, real64), &
                                   cmplx(5.0_real64, 6.0_real64, real64), &
                                   cmplx(7.0_real64, 8.0_real64, real64), &
                                   cmplx(9.0_real64, 10.0_real64, real64)]

    expected_header = [character(len=64) :: "IntColumn", "RealColumn", "CharColumn", "LogicalColumn", "ComplexColumn"]
    expected_metadata = reshape([1, 1, 2, 1, 3, 1, 4, 1, 5, 1], [2, 5])

    ! Read the CSV file
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    ! Check that reading was successful
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed")

    ! Test integer column getters
    call get_int_column_by_index(int_cols, metadata, 1, retrieved_int_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Integer getter by index should succeed")
    call assert_equal_array_int(retrieved_int_column, expected_int_cols(:, 1), 5, "Retrieved integer column should match")

    call get_int_column_by_name(int_cols, metadata, header, "IntColumn", retrieved_int_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Integer getter by name should succeed")
    call assert_equal_array_int(retrieved_int_column, expected_int_cols(:, 1), 5, "Retrieved integer column by name should match")

    ! Test real column getters
    call get_real_column_by_index(real_cols, metadata, 2, retrieved_real_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Real getter by index should succeed")
    call assert_equal_array_real(retrieved_real_column, expected_real_cols(:, 1), 5, 1.0e-10_real64, "Retrieved real column should match")

    call get_real_column_by_name(real_cols, metadata, header, "RealColumn", retrieved_real_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Real getter by name should succeed")
    call assert_equal_array_real(retrieved_real_column, expected_real_cols(:, 1), 5, 1.0e-10_real64, "Retrieved real column by name should match")

    ! Test character column getters
    call get_char_column_by_index(char_cols, metadata, 3, retrieved_char_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Character getter by index should succeed")
    call assert_equal_array_char(retrieved_char_column, expected_char_cols(:, 1), 64, 5, "Retrieved character column should match")

    call get_char_column_by_name(char_cols, metadata, header, "CharColumn", retrieved_char_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Character getter by name should succeed")
    call assert_equal_array_char(retrieved_char_column, expected_char_cols(:, 1), 64, 5, "Retrieved character column by name should match")

    ! Test logical column getters
    call get_logical_column_by_index(logical_cols, metadata, 4, retrieved_logical_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Logical getter by index should succeed")
    call assert_equal_array_int(merge(1, 0, retrieved_logical_column), merge(1, 0, expected_logical_cols(:, 1)), 5, "Retrieved logical column should match")
    
    call get_logical_column_by_name(logical_cols, metadata, header, "LogicalColumn", retrieved_logical_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Logical getter by name should succeed")
    call assert_equal_array_int(merge(1, 0, retrieved_logical_column), merge(1, 0, expected_logical_cols(:, 1)), 5, "Retrieved logical column should match")
   
    ! Test complex column getters
    call get_complex_column_by_index(complex_cols, metadata, 5, retrieved_complex_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Complex getter by index should succeed")
    call assert_equal_array_complex(retrieved_complex_column, expected_complex_cols(:, 1), 5, 1.0e-10_real64, "Retrieved complex column should match")

    call get_complex_column_by_name(complex_cols, metadata, header, "ComplexColumn", retrieved_complex_column, ierr)
    call assert_equal_int(ierr, ERR_OK, "Complex getter by name should succeed")
    call assert_equal_array_complex(retrieved_complex_column, expected_complex_cols(:, 1), 5, 1.0e-10_real64, "Retrieved complex column by name should match")

    ! Clean up temporary test file
    open(file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
  end subroutine test_getter_functions

  !> Test serialization and deserialization of CSV data
  subroutine test_serialization_deserialization()
    ! Define arrays for all 5 data types - mixed columns (5 rows)
    integer(int32) :: int_cols(5, 2), expected_int_cols(5, 2)
    real(real64) :: real_cols(5, 1), expected_real_cols(5, 1)
    character(len=64) :: char_cols(5, 2), expected_char_cols(5, 2)
    logical :: logical_cols(5, 1), expected_logical_cols(5, 1)
    complex(real64) :: complex_cols(5, 1), expected_complex_cols(5, 1)
    character(len=64) :: header(7), expected_header(7)
    integer(int32) :: metadata(2, 7), expected_metadata(2, 7)
    
    ! Arrays for deserialized data
    integer(int32) :: deserialized_int_cols(5, 2)
    real(real64) :: deserialized_real_cols(5, 1)
    character(len=64) :: deserialized_char_cols(5, 2)
    logical :: deserialized_logical_cols(5, 1)
    complex(real64) :: deserialized_complex_cols(5, 1)
    character(len=64) :: deserialized_header(7)
    integer(int32) :: deserialized_metadata(2, 7)
    
    integer(int32) :: ierr, i, file_unit, io_status
    character(len=*), parameter :: test_file = "test_serialize_temp.csv"
    character(len=*), parameter :: serialize_prefix = "testserial"
    
    ! Column types: int, int, real, char, char, logical, complex (7 columns total)
    integer(int32) :: column_types(7) = [1, 1, 2, 3, 3, 4, 5]
    
    ! Create temporary test CSV file with mixed data types
    open(newunit=file_unit, file=test_file, status='replace', action='write', iostat=io_status)
    if (is_err(io_status)) then
      print *, "Error creating test CSV file for serialization"
      stop 1
    end if
    
    ! Write CSV header
    write(file_unit, '(A)') "IntCol1,IntCol2,RealCol1,CharCol1,CharCol2,LogicalCol1,ComplexCol1"
    
    ! Write CSV data rows (5 rows with diverse data)
    write(file_unit, '(A)') "100,200,10.5,A,AA,T,(10.0,20.0)"
    write(file_unit, '(A)') "110,210,11.5,B,BB,F,(11.0,21.0)"
    write(file_unit, '(A)') "120,220,12.5,C,CC,T,(12.0,22.0)"
    write(file_unit, '(A)') "130,230,13.5,D,DD,F,(13.0,23.0)"
    write(file_unit, '(A)') "140,240,14.5,E,EE,T,(14.0,24.0)"

    close(file_unit)

    ! Set up expected data for verification
    expected_int_cols(:, 1) = [100, 110, 120, 130, 140]
    expected_int_cols(:, 2) = [200, 210, 220, 230, 240]
    expected_real_cols(:, 1) = [10.5_real64, 11.5_real64, 12.5_real64, 13.5_real64, 14.5_real64]
    expected_char_cols(:, 1) = [character(len=64) :: "A", "B", "C", "D", "E"]
    expected_char_cols(:, 2) = [character(len=64) :: "AA", "BB", "CC", "DD", "EE"]
    expected_logical_cols(:, 1) = [.true., .false., .true., .false., .true.]
    expected_complex_cols(:, 1) = [cmplx(10.0_real64, 20.0_real64, real64), &
                                   cmplx(11.0_real64, 21.0_real64, real64), &
                                   cmplx(12.0_real64, 22.0_real64, real64), &
                                   cmplx(13.0_real64, 23.0_real64, real64), &
                                   cmplx(14.0_real64, 24.0_real64, real64)]
    
    expected_header = [character(len=64) :: "IntCol1", "IntCol2", "RealCol1", "CharCol1", "CharCol2", "LogicalCol1", "ComplexCol1"]
    expected_metadata = reshape([1, 1, 1, 2, 2, 1, 3, 1, 3, 2, 4, 1, 5, 1], [2, 7])
    
    ! Read the CSV file into type-banded arrays
    call read_table(test_file, column_types, .true., int_cols, real_cols, char_cols, &
                    logical_cols, complex_cols, header, metadata, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "CSV reading should succeed for serialization test")
    
    ! Serialize the data to binary files
    call serialize_table(serialize_prefix, int_cols, real_cols, char_cols, &
                        logical_cols, complex_cols, header, metadata, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "Serialization should succeed")
    
    ! Deserialize the data back from binary files
    call deserialize_table(serialize_prefix, deserialized_int_cols, deserialized_real_cols, &
                          deserialized_char_cols, deserialized_logical_cols, &
                          deserialized_complex_cols, deserialized_header, deserialized_metadata, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "Deserialization should succeed")
    
    ! Verify that deserialized data matches original data
    call assert_equal_array_int(deserialized_int_cols, expected_int_cols, 10, "Deserialized integer columns should match original")
    call assert_equal_array_real(deserialized_real_cols, expected_real_cols, 5, 1.0e-10_real64, "Deserialized real columns should match original")
    call assert_equal_array_char(deserialized_char_cols, expected_char_cols, 64, 10, "Deserialized character columns should match original")
    call assert_equal_array_int(merge(1, 0, deserialized_logical_cols), merge(1, 0, expected_logical_cols), 5, "Deserialized logical columns should match original")
    call assert_equal_array_complex(deserialized_complex_cols, expected_complex_cols, 5, 1.0e-10_real64, "Deserialized complex columns should match original")
    call assert_equal_array_char(deserialized_header, expected_header, 64, 7, "Deserialized headers should match original")
    call assert_equal_array_int(deserialized_metadata, expected_metadata, 14, "Deserialized metadata should match original")
    
    ! Clean up temporary CSV file
    open(newunit=file_unit, file=test_file, status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    ! Clean up serialized binary files
    open(newunit=file_unit, file=serialize_prefix // "_int_cols.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_real_cols.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_char_cols.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_logical_cols.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_complex_cols.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_header.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if
    
    open(newunit=file_unit, file=serialize_prefix // "_metadata.dat", status='old', iostat=io_status)
    if (is_ok(io_status)) then
      close(file_unit, status='delete')
    end if

  end subroutine test_serialization_deserialization

end module mod_test_csv_file_reader