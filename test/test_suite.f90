! filepath: test/mod_test_suite.f90
!> Common definitions for all test suites.
module mod_test_suite
  implicit none
  private

  public :: test_interface
  public :: test_case
  public :: get_all_interface

  !> Common interface for all test procedures.
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  !> Type to hold information about a single test.
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

  !> Abstract interface for getting all tests in a suite.
  abstract interface
    function get_all_interface() result(all_tests)
      import :: test_case
      type(test_case), allocatable :: all_tests(:)
    end function get_all_interface
  end interface

end module mod_test_suite