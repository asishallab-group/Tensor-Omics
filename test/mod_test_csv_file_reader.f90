! filepath: test/mod_test_csv_file_reader.f90
!> Unit test suite for CSV file reader module.
module mod_test_csv_file_reader
  use asserts
  use tox_csv_file_reader
  use tox_errors, only: ERR_OK, ERR_EMPTY_INPUT, ERR_DIM_MISMATCH
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

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(1)
    all_tests(1) = test_case("test_sample_table", test_sample_table)
  end function get_all_tests

  !> Run all CSV file reader tests.
  subroutine run_all_tests_csv_file_reader()
    type(test_case) :: all_tests(1)
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
    type(test_case) :: all_tests(1)
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

  !> Test correct mapping between families and genes.
  subroutine test_sample_table()
    integer(int32) :: int_cols(2,1), i, j
    real(real64) :: real_cols(2,2)
    character(len=64) :: char_cols(2,1)
    logical :: logical_cols(2,1)
    complex(real64) :: complex_cols(2,5)
    character(len=64) :: header(1)
    integer(int32) :: metadata(2,5)
    integer(int32) :: ierr
    call read_table("test.csv", [1, 2, 3, 2, 4], .false., int_cols, real_cols, char_cols, &
                             logical_cols, complex_cols, header, metadata, ierr)
    
    print *, "Integer Columns:"
    do i = 1, 2
      do j = 1, 1
        print *, "Row ", i, ", Col ", j, ": ", int_cols(i, j)
      end do
    end do

    print *, "Real Columns:"
    do i = 1, 2
      do j = 1, 2
        print *, "Row ", i, ", Col ", j, ": ", real_cols(i, j)
      end do
    end do

    print *, "Character Columns:"
    do i = 1, 2
      do j = 1, 1
        print *, "Row ", i, ", Col ", j, ": ", char_cols(i, j)
      end do
    end do

    print *, "Logical Columns:"
    do i = 1, 2
      do j = 1, 1
        print *, "Row ", i, ", Col ", j, ": ", logical_cols(i, j)
      end do
    end do

    do i = 1, 2
      do j = 1, 5
        print *, "Metadata Row ", i, ", Col ", j, ": ", metadata(i, j)
      end do
    end do
  end subroutine test_sample_table

end module mod_test_csv_file_reader