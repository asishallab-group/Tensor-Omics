! filepath: test/test_suite.f90
!> Common definitions for all test suites.
module test_suite
    implicit none

    !> Common interface for all test procedures.
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    !> Type to hold information about a single test.
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    !> Abstract interface for getting all tests in a suite.
    abstract interface
        function get_all_interface() result(all_tests)
            import :: test_case
            type(test_case), allocatable :: all_tests(:)
        end function get_all_interface
    end interface

    !> Type to hold suite information.
    type :: suite_entry
        character(len=128) :: name
        procedure(get_all_interface), pointer, nopass :: get_all => null()
    end type suite_entry

    type(suite_entry), allocatable :: available_suites(:)

    !> Failed assertions in the test case currently running, and what the run has seen so far.
    !|
    !| A failed assertion used to end the process, so a run reported the first broken test and
    !| nothing about the rest -- which turns fixing a batch of related breakage into one
    !| rebuild per failure. It now fails the *case*: `run_test_case` checks this after the
    !| test returns, and the run keeps going. The remaining assertions of a broken test still
    !| execute, so a test that asserts a precondition and then relies on it can still crash the
    !| process -- but everything learned before that point is already on screen.
    integer, private :: failures_in_case = 0
    integer :: cases_run = 0
    integer :: cases_failed = 0
    integer :: suites_failed = 0

#ifndef NO_COLORS
    character(len=*), parameter :: COLOR_GREEN = achar(27) // "[38;5;154m"
    character(len=*), parameter :: COLOR_COPPER = achar(27) // "[38;5;214m"
    character(len=*), parameter :: COLOR_DARK_COPPER = achar(27) // "[38;5;208m"
    character(len=*), parameter :: COLOR_RED = achar(27) // "[38;5;196m"
    character(len=*), parameter :: COLOR_LIGHT_GRAY = achar(27) // "[38;5;252m"
    character(len=*), parameter :: COLOR_YELLOW = achar(27) // "[38;5;226m"
    character(len=*), parameter :: COLOR_CREAM = achar(27) // "[38;5;255m"
    character(len=*), parameter :: COLOR_ERROR = achar(27) // "[38;5;222m"
    character(len=*), parameter :: COLOR_RESET = achar(27) // "[0m"
#else
    character(len=*), parameter :: COLOR_GREEN = ""
    character(len=*), parameter :: COLOR_COPPER = ""
    character(len=*), parameter :: COLOR_DARK_COPPER = ""
    character(len=*), parameter :: COLOR_RED = ""
    character(len=*), parameter :: COLOR_LIGHT_GRAY = ""
    character(len=*), parameter :: COLOR_YELLOW = ""
    character(len=*), parameter :: COLOR_CREAM = ""
    character(len=*), parameter :: COLOR_ERROR = ""
    character(len=*), parameter :: COLOR_RESET = ""
#endif

    character(len=3) :: CHECK_MARK = char(int(z"E2")) // char(int(z"9C")) // char(int(z"93"))
    character(len=3) :: CROSS_MARK = char(int(z"E2")) // char(int(z"9C")) // char(int(z"97"))

contains

    !> Record a failed assertion against the test case currently running.
    !|
    !| Called by the assertion module, which cannot hold this state itself: it already `use`s
    !| this module for the colours, so the dependency only goes one way.
    subroutine record_assertion_failure()
        failures_in_case = failures_in_case + 1
    end subroutine record_assertion_failure

    !> Whether anything has failed so far.
    logical function any_failures()
        any_failures = cases_failed > 0
    end function any_failures

    !> Summarise the run and end the process with a status that reflects it.
    subroutine report_and_exit()
        character(len=:), allocatable :: counts

        if (cases_failed == 0) then
            print "(A)", COLOR_GREEN // "All " // COLOR_CREAM // itoa(cases_run) // &
                COLOR_GREEN // " test cases passed." // COLOR_RESET
            return
        end if

        counts = itoa(cases_failed) // " of " // itoa(cases_run) // " test cases failed"
        print "(A)", ""
        print "(A)", COLOR_RED // counts // COLOR_CREAM // &
            " (across " // itoa(suites_failed) // " suite(s))." // COLOR_RESET
        stop 1
    end subroutine report_and_exit

    !> An integer as a trimmed string, for the counts above.
    function itoa(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=32) :: buffer

        write (buffer, "(I0)") value
        text = trim(buffer)
    end function itoa

    !> Register all test suites here.
    subroutine initialize_suites()
        allocate (available_suites(0))
    end subroutine initialize_suites

    !> Add a new test suite to the registry.
    subroutine add_suite(name, get_all_proc)
        character(len=*), intent(in) :: name
        procedure(get_all_interface) :: get_all_proc
        type(suite_entry), allocatable :: temp_suites(:)
        integer :: n

        n = size(available_suites)
        allocate (temp_suites(n + 1))

        if (n > 0) then
            temp_suites(1:n) = available_suites(1:n)
        end if

        temp_suites(n + 1)%name = name
        temp_suites(n + 1)%get_all => get_all_proc

        call move_alloc(temp_suites, available_suites)
    end subroutine add_suite

    !> Run all tests in a given suite.
    subroutine run_all_tests(suite_name, all_tests)
        character(len=*), intent(in) :: suite_name
        type(test_case), intent(in) :: all_tests(:)
        integer :: i

        integer :: failed_before

        failed_before = cases_failed
        do i = 1, size(all_tests)
            call run_test_case(all_tests(i))
        end do

        if (cases_failed == failed_before) then
            print "(A)", COLOR_CREAM // "All '" // COLOR_LIGHT_GRAY // trim(suite_name) // COLOR_CREAM // "' tests passed successfully." // COLOR_RESET
        else
            suites_failed = suites_failed + 1
            print "(A)", COLOR_RED // itoa(cases_failed - failed_before) // COLOR_CREAM // " test(s) in '" // &
                COLOR_LIGHT_GRAY // trim(suite_name) // COLOR_CREAM // "' failed." // COLOR_RESET
        end if
    end subroutine run_all_tests

    !> Run test for specific test case
    subroutine run_test_case(case)
        type(test_case), intent(in) :: case

        failures_in_case = 0
        cases_run = cases_run + 1

        call case%test_proc()

        if (failures_in_case == 0) then
            print "(A)", COLOR_GREEN // CHECK_MARK // COLOR_COPPER // " " // trim(case%name) // COLOR_GREEN // " passed" // COLOR_CREAM // "." // COLOR_RESET
        else
            cases_failed = cases_failed + 1
            print "(A)", COLOR_RED // CROSS_MARK // COLOR_COPPER // " " // trim(case%name) // COLOR_RED // " failed" // &
                COLOR_CREAM // " (" // itoa(failures_in_case) // " assertion(s))." // COLOR_RESET
        end if
    end subroutine run_test_case

    !> Run selected tests by name from a given suite.
    subroutine run_named_tests(test_names, all_tests)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), intent(in) :: all_tests(:)
        integer :: i, j
        logical :: found, some_not_found

        some_not_found = .false.
        do i = 1, size(test_names)
            found = .false.

            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call run_test_case(all_tests(j))
                    found = .true.
                    exit
                end if
            end do

            if (.not. found) then
                print "(A)", COLOR_RED // "Unknown test" // COLOR_CREAM // ": " // COLOR_COPPER // trim(test_names(i)) // COLOR_RESET
                some_not_found = .true.
            end if
        end do

        ! It is an error if a test case doesn't exist
        if (some_not_found) then
            stop 1
        end if
    end subroutine run_named_tests

    !> Run every registered suite.
    subroutine run_all_suites()
        integer :: i
        type(test_case), allocatable :: all_tests(:)

        do i = 1, size(available_suites)
            print "(A)", COLOR_CREAM // "Running suite: '" // COLOR_LIGHT_GRAY // trim(available_suites(i)%name) // COLOR_CREAM // "'" // COLOR_RESET
            all_tests = available_suites(i)%get_all()
            call run_all_tests(trim(available_suites(i)%name), all_tests)
        end do
    end subroutine run_all_suites

    !> Run one whole suite by name.
    subroutine run_suite_all(requested_suite)
        character(len=*), intent(in) :: requested_suite
        integer :: i
        type(test_case), allocatable :: all_tests(:)

        do i = 1, size(available_suites)
            if (trim(available_suites(i)%name) == trim(requested_suite)) then
                print "(A)", COLOR_CREAM // "Running suite: '" // COLOR_LIGHT_GRAY // trim(available_suites(i)%name) // COLOR_CREAM // "'" // COLOR_RESET
                all_tests = available_suites(i)%get_all()
                call run_all_tests(trim(requested_suite), all_tests)
                return
            end if
        end do

        print "(A)", COLOR_RED // "Unknown test suite" // COLOR_CREAM // ": " // COLOR_DARK_COPPER // trim(requested_suite) // COLOR_RESET
        stop 1
    end subroutine run_suite_all

    !> Run selected tests from one suite.
    subroutine run_suite_named(requested_suite, test_list)
        character(len=*), intent(in) :: requested_suite, test_list
        integer :: i
        type(test_case), allocatable :: all_tests(:)
        character(len=128), allocatable :: test_names(:)

        do i = 1, size(available_suites)
            if (trim(available_suites(i)%name) == trim(requested_suite)) then
                all_tests = available_suites(i)%get_all()
                call split_test_list(test_list, test_names)
                call run_named_tests(test_names, all_tests)
                return
            end if
        end do

        print "(A)", COLOR_RED // "Unknown test suite" // COLOR_CREAM // ": " // COLOR_DARK_COPPER // trim(requested_suite) // COLOR_RESET
        stop 1
    end subroutine run_suite_named

    !> Split a comma-separated list of test names into an array.
    subroutine split_test_list(test_list, test_names)
        character(len=*), intent(in) :: test_list
        character(len=128), allocatable, intent(out) :: test_names(:)

        integer :: i, n_items, lenstr
        integer :: start_pos, comma_pos, item_idx

        lenstr = len_trim(test_list)

        if (lenstr == 0) then
            allocate (test_names(0))
            return
        end if

        n_items = 1
        do i = 1, lenstr
            if (test_list(i:i) == ',') n_items = n_items + 1
        end do

        allocate (test_names(n_items))
        test_names = ""

        start_pos = 1
        item_idx = 1

        do
            comma_pos = index(test_list(start_pos:), ',')

            if (comma_pos == 0) then
                test_names(item_idx) = trim(adjustl(test_list(start_pos:lenstr)))
                exit
            else
                test_names(item_idx) = trim(adjustl(test_list(start_pos:start_pos + comma_pos - 2)))
                start_pos = start_pos + comma_pos
                item_idx = item_idx + 1
            end if
        end do
    end subroutine split_test_list

end module test_suite
